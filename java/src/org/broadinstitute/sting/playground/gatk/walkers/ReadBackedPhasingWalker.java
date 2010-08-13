/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriterImpl;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.*;
import java.util.*;


/**
 * Walks along all loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using downstream reads).
 * Use '-BTI variant' to only stop at positions in the VCF file bound to 'variant'.
 */
@Requires(value = {}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))
public class ReadBackedPhasingWalker extends LocusWalker<PhasingStatsAndOutput, PhasingStats> {

    @Argument(fullName = "cacheWindowSize", shortName = "cacheWindow", doc = "The window size (in bases) to cache variant sites and their reads; [default:20000]", required = false)
    protected Integer cacheWindow = 20000;

    @Argument(fullName = "phasedVCFFile", shortName = "phasedVCF", doc = "The name of the phased VCF file output", required = true)
    protected String phasedVCFFile = null;

    private VCFWriter writer = null;

    private LinkedList<VariantAndReads> siteQueue = null;
    private VariantAndReads prevVr = null; // the VC emitted after phasing, and the alignment bases at the position emitted

    private static double SMALL_THRESH = 1e-6;
    private static int MAX_NUM_PHASE_SITES = 20; // 2^20 == 10^6 biallelic haplotypes

    private void initializeVcfWriter(VariantContext vc) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.add(new VCFFormatHeaderLine("PQ", 1, VCFHeaderLineType.Float, "Read-backed phasing quality score"));

        writer = new VCFWriterImpl(new File(phasedVCFFile));
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(vc.getSampleNames())));
    }

    public void initialize() {
        siteQueue = new LinkedList<VariantAndReads>();
        prevVr = new VariantAndReads(null, null, true);
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public PhasingStats reduceInit() {
        return new PhasingStats();
    }

    /**
     * For each site of interest, cache the current site and then use the cache to phase all upstream sites
     * for which "sufficient" information has already been observed.
     *
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public PhasingStatsAndOutput map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        PhasingStats phaseStats = new PhasingStats();

        LinkedList<String> rodNames = new LinkedList<String>();
        rodNames.add("variant");
        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = false;
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            boolean processVariant = true;
            if (!vc.isSNP() || !vc.isBiallelic()) {
                processVariant = false;
            }
            VariantAndReads vr = new VariantAndReads(vc, context, processVariant);
            siteQueue.add(vr);

            int numReads = 0;
            if (context.hasBasePileup()) {
                numReads = context.getBasePileup().size();
            }
            else if (context.hasExtendedEventPileup()) {
                numReads = context.getExtendedEventPileup().size();
            }
            PhasingStats addInPhaseStats = new PhasingStats(numReads, 1);
            phaseStats.addIn(addInPhaseStats);
        }
        List<VariantContext> phasedList = processQueue(ref.getLocus(), phaseStats);

        return new PhasingStatsAndOutput(phaseStats, phasedList);
    }

    private List<VariantContext> processQueue(GenomeLoc loc, PhasingStats phaseStats) {
        List<VariantContext> vcList = new LinkedList<VariantContext>();

        while (!siteQueue.isEmpty()) {
            if (loc != null) {
                VariantContext vc = siteQueue.peek().variant;
                if (isInWindowRange(loc, VariantContextUtils.getLocation(vc))) {
                    // loc is still not far enough ahead of vc (since we ASSUME that the VCF is ordered by <contig,locus>)
                    break;
                }
                // Already saw all variant positions within cacheWindow distance ahead of vc (on its contig)
            }
            VariantContext phasedVc = finalizePhasingAndRemove(phaseStats);
            vcList.add(phasedVc);
        }
        return vcList;
    }

    /* Finalize phasing of vc (head of siteQueue) using all VariantContext objects in the siteQueue that are
        within cacheWindow distance ahead of vc (on its contig).
        ASSUMES:
        1. siteQueue is NOT empty.
        2. All VariantContexts in siteQueue are in positions downstream of vc (head of queue).
     */

    private VariantContext finalizePhasingAndRemove(PhasingStats phaseStats) {
        VariantAndReads vr = siteQueue.remove(); // remove vr from head of queue
        VariantContext vc = vr.variant;
        if (!vr.processVariant)
            return vc; // return vc as is

        boolean hasPreviousSite = previousIsRelevantTo(vc);

        logger.debug("Will phase vc = " + VariantContextUtils.getLocation(vc));

        LinkedList<VariantAndReads> windowVaList = new LinkedList<VariantAndReads>();
        if (hasPreviousSite) {
            windowVaList.add(prevVr); // need to add one position for phasing context
            windowVaList.add(vr); // add position to be phased
            for (VariantAndReads nextVr : siteQueue) {
                if (!isInWindowRange(vc, nextVr.variant)) //nextVr too far ahead of the range used for phasing vc
                    break;
                if (nextVr.processVariant) // include in the phasing computation
                    windowVaList.add(nextVr);
            }

            if (logger.isDebugEnabled()) {
                ListIterator<VariantAndReads> windowVcIt = windowVaList.listIterator();
                while (windowVcIt.hasNext()) {
                    VariantContext phaseInfoVc = windowVcIt.next().variant;
                    logger.debug("Using phaseInfoVc = " + VariantContextUtils.getLocation(phaseInfoVc));
                }
            }
        }

        logger.debug("");

        Map<String, Genotype> sampGenotypes = vc.getGenotypes();
        VariantContext prevVc = prevVr.variant;
        Map<String, Genotype> phasedGtMap = new TreeMap<String, Genotype>();

        // Perform per-sample phasing:
        TreeMap<String, PhaseCounts> samplePhaseStats = new TreeMap<String, PhaseCounts>();
        for (Map.Entry<String, Genotype> sampGtEntry : sampGenotypes.entrySet()) {
            logger.debug("sample = " + sampGtEntry.getKey());

            boolean genotypesArePhased = true; // phase by default

            String samp = sampGtEntry.getKey();
            Genotype gt = sampGtEntry.getValue();
            Biallele biall = new Biallele(gt);
            HashMap<String, Object> gtAttribs = new HashMap<String, Object>(gt.getAttributes());

            if (hasPreviousSite && gt.isHet() && prevVc.getGenotype(samp).isHet()) { //otherwise, can trivially phase
                logger.debug("NON-TRIVIALLY CARE about TOP vs. BOTTOM for: ");
                logger.debug("\n" + biall);

                LinkedList<VariantAndReads> sampleWindowVaList = new LinkedList<VariantAndReads>();
                for (VariantAndReads phaseInfoVr : windowVaList) {
                    VariantContext phaseInfoVc = phaseInfoVr.variant;
                    Genotype phaseInfoGt = phaseInfoVc.getGenotype(samp);
                    if (phaseInfoGt.isHet()) { // otherwise, of no value to phasing
                        sampleWindowVaList.add(phaseInfoVr);
                        logger.debug("STARTING TO PHASE USING POS = " + VariantContextUtils.getLocation(phaseInfoVc));
                    }
                }
                if (sampleWindowVaList.size() > MAX_NUM_PHASE_SITES)
                    logger.warn("Trying to phase within a window of " + cacheWindow + " bases yields " + sampleWindowVaList.size() + " heterozygous sites to phase -- EXPECT DELAYS!");

                PhasingTable sampleHaps = new PhasingTable();

                // Initialize phasing table with appropriate entries:
                //
                // 1. THIS IMPLEMENTATION IS INEFFICIENT SINCE IT DOES NOT PREALLOCATE
                // THE ArrayList USED, BUT RATHER APPENDS TO IT EACH TIME.
                //
                // 2. THIS IMPLEMENTATION WILL FAIL WHEN NOT DEALING WITH SNP Alleles, SINCE THEN THE Allele.getBases()
                // FUNCTION WILL RETURN VARIABLE-LENGTH Byte ARRAYS.  IN THAT CASE, BaseArray WILL NEED TO BE CONVERTED TO
                // AN ArrayList OF Allele [OR SIMILAR OBJECT]
                for (VariantAndReads phaseInfoVr : sampleWindowVaList) {
                    VariantContext phaseInfoVc = phaseInfoVr.variant;
                    Genotype phaseInfoGt = phaseInfoVc.getGenotype(samp);

                    if (sampleHaps.isEmpty()) {
                        for (Allele sampAll : phaseInfoGt.getAlleles()) {
                            sampleHaps.addEntry(new Haplotype(sampAll.getBases()));
                        }
                    }
                    else {
                        PhasingTable oldHaps = sampleHaps;
                        Iterator<PhasingTable.PhasingTableEntry> oldHapIt = oldHaps.iterator();
                        sampleHaps = new PhasingTable();
                        while (oldHapIt.hasNext()) {
                            PhasingTable.PhasingTableEntry pte = oldHapIt.next();
                            Haplotype oldHap = pte.getHaplotype();
                            for (Allele sampAll : phaseInfoGt.getAlleles()) {
                                ArrayList<Byte> bases = oldHap.cloneBaseArrayList();
                                for (byte b : sampAll.getBases()) { // LENGTH NOT PRE-DEFINED FOR NON-SNPs (MNP or INDELS!!)
                                    bases.add(b); // INEFFICIENT!
                                }
                                Haplotype newHap = new Haplotype(BaseArray.getBasesPrimitiveNoNulls(bases));
                                sampleHaps.addEntry(newHap);
                            }
                        }
                    }
                }

                // Assemble the "sub-reads" from the heterozygous positions for this sample:
                LinkedList<ReadBasesAtPosition> allPositions = new LinkedList<ReadBasesAtPosition>();
                for (VariantAndReads phaseInfoVr : sampleWindowVaList) {
                    ReadBasesAtPosition readBases = phaseInfoVr.sampleReadBases.get(samp);
                    allPositions.add(readBases);
                }
                HashMap<String, Read> allReads = convertReadBasesAtPositionToReads(allPositions);
                logger.debug("Number of reads at sites: " + allReads.size());

                // Update the phasing table based on each of the sub-reads for this sample:
                int numUsableReads = 0;
                for (Map.Entry<String, Read> nameToReads : allReads.entrySet()) {
                    Read rd = nameToReads.getValue();
                    if (rd.numNonNulls() <= 1) {// can't possibly provide any phasing information
                        continue;
                    }
                    if (false)
                        logger.debug("rd = " + rd + "\tname = " + nameToReads.getKey());

                    LinkedList<PhasingTable.PhasingTableEntry> compatHaps = new LinkedList<PhasingTable.PhasingTableEntry>();
                    Iterator<PhasingTable.PhasingTableEntry> hapIt = sampleHaps.iterator();
                    while (hapIt.hasNext()) {
                        PhasingTable.PhasingTableEntry pte = hapIt.next();
                        if (rd.isCompatible(pte.getHaplotype()))
                            compatHaps.add(pte);
                    }

                    if (!compatHaps.isEmpty()) { // otherwise, nothing to do
                        numUsableReads++;
                        double addVal = rd.matchScore() / compatHaps.size(); // don't overcount, so divide up the score evenly
                        for (PhasingTable.PhasingTableEntry pte : compatHaps) {
                            pte.addScore(addVal);

                            if (false) {
                                if (addVal > SMALL_THRESH) {
                                    logger.debug("score(" + rd + "," + pte.getHaplotype() + ") = " + addVal);
                                }
                            }
                        }
                    }
                }
                logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");
                logger.debug("numUsableReads = " + numUsableReads);

                /* Map a phase and its "complement" to a single representative phase, but marginalized to the first 2 positions
                  [i.e., the previous position and the current position]:
                 */
                ComplementAndMarginalizeHaplotypeMapper cmhm = new ComplementAndMarginalizeHaplotypeMapper(sampleWindowVaList, samp, 2);
                sampleHaps = sampleHaps.mapHaplotypes(cmhm);

                logger.debug("\nPhasing table [AFTER MAPPING]:\n" + sampleHaps + "\n");

                // Determine the phase at this position:
                sampleHaps.normalizeScores();
                PhasingTable.PhasingTableEntry maxEntry = sampleHaps.maxEntry();
                double score = maxEntry.getScore();
                genotypesArePhased = (score > SMALL_THRESH);

                if (genotypesArePhased) {
                    Biallele prevBiall = new Biallele(prevVc.getGenotype(samp));
                    ensurePhasing(maxEntry.getHaplotype(), biall, prevBiall);
                    gtAttribs.put("PQ", new Float(score));

                    logger.debug("CHOSE hap:\t" + maxEntry.getHaplotype() + "\tscore:\t" + score);
                    logger.debug("PHASED:\n" + biall + "\n\n");
                }

                PhaseCounts sampPhaseCounts = samplePhaseStats.get(samp);
                if (sampPhaseCounts == null) {
                    sampPhaseCounts = new PhaseCounts();
                    samplePhaseStats.put(samp, sampPhaseCounts);
                }
                sampPhaseCounts.numTestedSites++;
                sampPhaseCounts.numPhased += (genotypesArePhased ? 1 : 0);
            }

            List<Allele> phasedAll = biall.getAllelesAsList();
            Genotype phasedGt = new Genotype(gt.getSampleName(), phasedAll, gt.getNegLog10PError(), gt.getFilters(), gtAttribs, genotypesArePhased);
            phasedGtMap.put(samp, phasedGt);
        }

        VariantContext phasedVc = new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), phasedGtMap, vc.getNegLog10PError(), vc.getFilters(), vc.getAttributes());
        prevVr.variant = phasedVc;
        prevVr.sampleReadBases = vr.sampleReadBases;

        phaseStats.addIn(new PhasingStats(samplePhaseStats));

        return phasedVc;
    }

    /*
        Ensure that curBiall is phased relative to prevBiall as specified by hap.
     */

    public static void ensurePhasing(Haplotype hap, Biallele curBiall, Biallele prevBiall) {
        if (hap.size() < 2)
            throw new StingException("LOGICAL ERROR: Only considering haplotypes of length > 2!");

        byte prevBase = hap.getBase(0); // The 1st base in the haplotype
        byte curBase = hap.getBase(1);  // The 2nd base in the haplotype

        boolean chosePrevTopChrom = prevBiall.matchesTopBase(prevBase);
        boolean choseCurTopChrom = curBiall.matchesTopBase(curBase);
        if (chosePrevTopChrom != choseCurTopChrom)
            curBiall.swapAlleles();
    }

    private boolean previousIsRelevantTo(VariantContext vc) {
        VariantContext prevVc = prevVr.variant;
        return (prevVc != null && VariantContextUtils.getLocation(prevVc).onSameContig(VariantContextUtils.getLocation(vc)));
    }

    private boolean isInWindowRange(VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(vc2);

        return isInWindowRange(loc1, loc2);
    }

    private boolean isInWindowRange(GenomeLoc loc1, GenomeLoc loc2) {
        return (loc1.onSameContig(loc2) && loc1.distance(loc2) <= cacheWindow);
    }

    private void writeVCF(VariantContext vc) {
        if (writer == null)
            initializeVcfWriter(vc);

        byte refBase;
        if (!vc.isIndel()) {
            Allele varAllele = vc.getReference();
            refBase = varAllele.getBases()[0];
        }
        else {
            refBase = vc.getReferenceBaseForIndel();
        }

        writer.add(vc, refBase);
    }

    public PhasingStats reduce(PhasingStatsAndOutput statsAndList, PhasingStats stats) {
        if (statsAndList != null) {
            writeVarContIter(statsAndList.output.iterator());
            stats.addIn(statsAndList.ps);
        }
        return stats;
    }

    /**
     * Phase anything left in the cached siteQueue, and report the number of reads and VariantContexts processed.
     *
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(PhasingStats result) {
        List<VariantContext> finalList = processQueue(null, result);
        writeVarContIter(finalList.iterator());
        if (writer != null)
            writer.close();

        out.println("Number of reads observed: " + result.getNumReads());
        out.println("Number of variant sites observed: " + result.getNumVarSites());
        out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));

        out.println("\n-- Phasing summary --");
        for (Map.Entry<String, PhaseCounts> sampPhaseCountEntry : result.getPhaseCounts()) {
            PhaseCounts pc = sampPhaseCountEntry.getValue();
            out.println("Sample: " + sampPhaseCountEntry.getKey() + "\tNumber of tested sites: " + pc.numTestedSites + "\tNumber of phased sites: " + pc.numPhased);
        }
        out.println("");
    }

    protected void writeVarContIter(Iterator<VariantContext> varContIter) {
        while (varContIter.hasNext()) {
            VariantContext vc = varContIter.next();
            writeVCF(vc);
        }
    }

    protected static HashMap<String, Read> convertReadBasesAtPositionToReads(Collection<ReadBasesAtPosition> basesAtPositions) {
        HashMap<String, Read> reads = new HashMap<String, Read>();

        int index = 0;
        for (ReadBasesAtPosition rbp : basesAtPositions) {
            Iterator<ReadBasesAtPosition.ReadBase> readBaseIt = rbp.iterator();
            while (readBaseIt.hasNext()) {
                ReadBasesAtPosition.ReadBase rb = readBaseIt.next();
                String readName = rb.readName;
                byte base = rb.base;

                Read rd = reads.get(readName);
                if (rd == null) {
                    rd = new Read(basesAtPositions.size());
                    reads.put(readName, rd);
                }
                rd.updateBase(index, base);
            }
            index++;
        }

        return reads;
    }


    /*
       Inner classes:
     */

    private static class Biallele {
        public Allele top;
        public Allele bottom;

        public Biallele(Genotype gt) {
            if (gt.getPloidy() != 2)
                throw new StingException("Doesn't support phasing for ploidy that is not 2!");

            this.top = gt.getAllele(0);
            this.bottom = gt.getAllele(1);
        }

        public void swapAlleles() {
            Allele tmp = top;
            top = bottom;
            bottom = tmp;
        }

        public List<Allele> getAllelesAsList() {
            List<Allele> allList = new ArrayList<Allele>(2);
            allList.add(0, top);
            allList.add(1, bottom);
            return allList;
        }

        public byte getTopBase() {
            byte[] topBases = top.getBases();
            if (topBases.length != 1)
                throw new StingException("LOGICAL ERROR: should not process non-SNP sites!");

            return topBases[0];
        }

        public byte getBottomBase() {
            byte[] bottomBases = bottom.getBases();
            if (bottomBases.length != 1)
                throw new StingException("LOGICAL ERROR: should not process non-SNP sites!");

            return bottomBases[0];
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Top:\t" + top.getBaseString() + "\n");
            sb.append("Bot:\t" + bottom.getBaseString() + "\n");
            return sb.toString();
        }

        public boolean matchesTopBase(byte base) {
            boolean matchesTop;
            if (base == getTopBase())
                matchesTop = true;
            else if (base == getBottomBase())
                matchesTop = false;
            else
                throw new StingException("LOGICAL ERROR: base MUST match either TOP or BOTTOM!");

            return matchesTop;
        }

        public byte getOtherBase(byte base) {
            byte topBase = getTopBase();
            byte botBase = getBottomBase();

            if (base == topBase)
                return botBase;
            else if (base == botBase)
                return topBase;
            else
                throw new StingException("LOGICAL ERROR: base MUST match either TOP or BOTTOM!");
        }
    }

    private static class VariantAndReads {
        public VariantContext variant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;
        public boolean processVariant;

        public VariantAndReads(VariantContext variant, AlignmentContext alignment, boolean processVariant) {
            this.variant = variant;
            this.sampleReadBases = new HashMap<String, ReadBasesAtPosition>();
            this.processVariant = processVariant;

            if (alignment != null) {
                ReadBackedPileup pileup = null;
                if (alignment.hasBasePileup()) {
                    pileup = alignment.getBasePileup();
                }
                else if (alignment.hasExtendedEventPileup()) {
                    pileup = alignment.getExtendedEventPileup();
                }
                if (pileup != null) {
                    for (String samp : pileup.getSamples()) {
                        ReadBackedPileup samplePileup = pileup.getPileupForSample(samp);
                        ReadBasesAtPosition readBases = new ReadBasesAtPosition();
                        for (PileupElement p : samplePileup) {
                            if (!p.isDeletion())
                                readBases.putReadBase(p.getRead().getReadName(), p.getBase());
                        }
                        sampleReadBases.put(samp, readBases);
                    }
                }
            }
        }
    }

    private static class ReadBasesAtPosition {
        // list of: <read name, base>
        private LinkedList<ReadBase> bases;

        public ReadBasesAtPosition() {
            this.bases = new LinkedList<ReadBase>();
        }

        public void putReadBase(String readName, byte b) {
            bases.add(new ReadBase(readName, b));
        }

        public Iterator<ReadBase> iterator() {
            return bases.iterator();
        }

        private static class ReadBase {
            public String readName;
            public byte base;

            public ReadBase(String readName, byte base) {
                this.readName = readName;
                this.base = base;
            }
        }
    }

    private static abstract class HaplotypeMapper {
        abstract Haplotype map(Haplotype hap);
    }

    private static class ComplementAndMarginalizeHaplotypeMapper extends HaplotypeMapper {
        private List<VariantAndReads> vaList;
        private String sample;
        private int marginalizeLength;

        public ComplementAndMarginalizeHaplotypeMapper(List<VariantAndReads> vaList, String sample, int marginalizeLength) {
            this.vaList = vaList;
            this.sample = sample;
            this.marginalizeLength = marginalizeLength;
        }

        public Haplotype map(Haplotype hap) {
            if (hap.size() != vaList.size())
                throw new StingException("INTERNAL ERROR: hap.size() != vaList.size()");

            Biallele firstPosBiallele = new Biallele(vaList.get(0).variant.getGenotype(sample));
            if (firstPosBiallele.matchesTopBase(hap.getBase(0))) {
                /* hap already matches the representative haplotype [arbitrarily defined to be
                   the one with the top base in the VariantContext at the 1st position]:
                 */
                return hap.subHaplotype(0, marginalizeLength); // only want first marginalizeLength positions
            }

            if (false)
                logger.debug("hap = " + hap);

            // Take the other base at EACH position of the Haplotype:
            byte[] complementBases = new byte[Math.min(hap.size(), marginalizeLength)];
            int index = 0;
            for (VariantAndReads vr : vaList) {
                VariantContext vc = vr.variant;
                Biallele biall = new Biallele(vc.getGenotype(sample));

                if (false)
                    logger.debug("biall =\n" + biall);

                complementBases[index] = biall.getOtherBase(hap.getBase(index));
                if (++index == marginalizeLength) // only want first marginalizeLength positions
                    break;
            }
            return new Haplotype(complementBases);
        }
    }

    private static class PhasingTable {
        private LinkedList<PhasingTableEntry> table;

        public PhasingTable() {
            this.table = new LinkedList<PhasingTableEntry>();
        }

        public PhasingTableEntry addEntry(Haplotype haplotype) {
            PhasingTableEntry pte = new PhasingTableEntry(haplotype);
            table.add(pte);
            return pte;
        }

        public Iterator<PhasingTableEntry> iterator() {
            return table.iterator();
        }

        public boolean isEmpty() {
            return table.isEmpty();
        }

        public PhasingTableEntry maxEntry() {
            if (table.isEmpty())
                return null;

            PhasingTableEntry maxPte = null;
            for (PhasingTableEntry pte : table) {
                if (maxPte == null || pte.getScore() > maxPte.getScore()) {
                    maxPte = pte;
                }
            }
            return maxPte.clone();
        }

        // Assumes that scores are NON-NEGATIVE:

        public void normalizeScores() {
            double normalizeBy = 0.0;
            for (PhasingTableEntry pte : table) {
                normalizeBy += pte.getScore();
            }
            logger.debug("normalizeBy = " + normalizeBy);

            if (normalizeBy > SMALL_THRESH) { // otherwise, will have precision problems
                for (PhasingTableEntry pte : table) {
                    pte.setScore(pte.getScore() / normalizeBy);
                }
            }
        }

        public PhasingTable mapHaplotypes(HaplotypeMapper hm) {
            class Score {
                private double d;

                Score(double d) {
                    this.d = d;
                }

                Score addValue(double v) {
                    d += v;
                    return this;
                }

                double value() {
                    return d;
                }
            }
            TreeMap<Haplotype, Score> hapMap = new TreeMap<Haplotype, Score>();

            Iterator<PhasingTableEntry> entryIt = iterator();
            while (entryIt.hasNext()) {
                PhasingTableEntry pte = entryIt.next();
                Haplotype rep = hm.map(pte.getHaplotype());
                if (false)
                    logger.debug("MAPPED: " + pte.getHaplotype() + " -> " + rep);

                Score score = hapMap.get(rep);
                if (score == null) {
                    score = new Score(0.0);
                    hapMap.put(rep, score);
                }
                score.addValue(pte.getScore());
            }

            PhasingTable combo = new PhasingTable();
            for (Map.Entry<Haplotype, Score> hapScore : hapMap.entrySet()) {
                combo.addEntry(hapScore.getKey()).addScore(hapScore.getValue().value());
            }
            return combo;
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("-------------------\n");
            Iterator<PhasingTable.PhasingTableEntry> hapIt = iterator();
            while (hapIt.hasNext()) {
                PhasingTable.PhasingTableEntry pte = hapIt.next();
                sb.append("Haplotype:\t" + pte.getHaplotype() + "\tScore:\t" + pte.getScore() + "\n");
            }
            sb.append("-------------------\n");
            return sb.toString();
        }

        public static class PhasingTableEntry implements Comparable<PhasingTableEntry>, Cloneable {
            private Haplotype haplotype;
            private double score;

            public PhasingTableEntry(Haplotype haplotype) {
                this.haplotype = haplotype;
                this.score = 0.0;
            }

            public PhasingTableEntry(PhasingTableEntry other) {
                this.haplotype = other.haplotype.clone();
                this.score = other.score;
            }

            public PhasingTableEntry clone() {
                try {
                    super.clone();
                } catch (CloneNotSupportedException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
                return new PhasingTableEntry(this);
            }

            public Haplotype getHaplotype() {
                return haplotype;
            }

            public double getScore() {
                return score;
            }

            public double addScore(double addVal) {
                score += addVal;
                return score;
            }

            private double setScore(double newVal) {
                score = newVal;
                return score;
            }

            public int compareTo(PhasingTableEntry that) {
                return new Double(Math.signum(this.score - that.score)).intValue();
            }
        }
    }
}


abstract class BaseArray implements Comparable<BaseArray> {
    protected ArrayList<Byte> bases;

    public BaseArray(byte[] bases) {
        this.bases = new ArrayList<Byte>(bases.length);
        for (int i = 0; i < bases.length; i++) {
            this.bases.add(i, bases[i]);
        }
    }

    public BaseArray(Byte[] bases) {
        this.bases = new ArrayList<Byte>(bases.length);
        for (int i = 0; i < bases.length; i++) {
            this.bases.add(i, bases[i]);
        }
    }

    public BaseArray(int length) {
        this(newNullArray(length));
    }

    static Byte[] newNullArray(int length) {
        Byte[] bArr = new Byte[length];
        Arrays.fill(bArr, null);
        return bArr;
    }

    public void updateBase(int index, Byte base) {
        bases.set(index, base);
    }

    public Byte getBase(int index) {
        return bases.get(index);
    }

    public int size() {
        return bases.size();
    }

    public static Byte[] getBases(List<Byte> baseList) {
        return baseList.toArray(new Byte[baseList.size()]);
    }

    // Will thow NullPointerException if baseList contains Byte == null:

    public static byte[] getBasesPrimitiveNoNulls(List<Byte> baseList) {
        int sz = baseList.size();
        byte[] b = new byte[sz];
        for (int i = 0; i < sz; i++) {
            b[i] = baseList.get(i);
        }
        return b;
    }

    public ArrayList<Byte> cloneBaseArrayList() {
        return new ArrayList<Byte>(bases);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(bases.size());
        for (Byte b : bases) {
            sb.append(b != null ? (char) b.byteValue() : "_");
        }
        return sb.toString();
    }

    public int compareTo(BaseArray that) {
        int sz = this.bases.size();
        if (sz != that.bases.size())
            return (sz - that.bases.size());

        for (int i = 0; i < sz; i++) {
            Byte thisBase = this.getBase(i);
            Byte thatBase = that.getBase(i);
            if (thisBase == null || thatBase == null) {
                if (thisBase == null && thatBase != null) {
                    return -1;
                }
                else if (thisBase != null && thatBase == null) {
                    return 1;
                }
            }
            else if (!thisBase.equals(thatBase)) {
                return thisBase - thatBase;
            }
        }
        return 0;
    }
}

class Haplotype extends BaseArray implements Cloneable {
    public Haplotype(byte[] bases) {
        super(bases);
    }

    private Haplotype(Byte[] bases) {
        super(bases);
    }

    private Haplotype(int length) {
        super(length);
    }

    public Haplotype(Haplotype other) {
        this(getBases(other.bases));
    }

    public void updateBase(int index, Byte base) {
        if (base == null) {
            throw new StingException("Internal error: should NOT put null for a missing Haplotype base!");
        }
        super.updateBase(index, base);
    }

    public Haplotype clone() {
        try {
            super.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return new Haplotype(this);
    }

    // Returns a new Haplotype containing the portion of this Haplotype between the specified fromIndex, inclusive, and toIndex, exclusive.

    public Haplotype subHaplotype(int fromIndex, int toIndex) {
        return new Haplotype(BaseArray.getBasesPrimitiveNoNulls(cloneBaseArrayList().subList(fromIndex, toIndex)));
    }
}

class Read extends BaseArray {
    //
    // ADD IN SOME DATA MEMBERS [OR INPUT TO CONSTRUCTORS] WITH READ QUALITY (MAPPING QUALITY) AND ARRAY OF BASE QUALITIES...
    // THESE WILL BE USED IN matchScore()
    //

    public Read(byte[] bases) {
        super(bases);
    }

    public Read(Byte[] bases) {
        super(bases);
    }

    public Read(int length) {
        super(length);
    }

    public int numNonNulls() {
        int num = 0;
        for (int i = 0; i < bases.size(); i++) {
            if (getBase(i) != null)
                num++;
        }
        return num;
    }

    //

    public double matchScore() {
        return 1.0;
    }
    //

    /* Checks if the two BaseArrays are consistent where bases are not null.
    */

    public boolean isCompatible(Haplotype hap) {
        int sz = this.bases.size();
        if (sz != hap.bases.size())
            throw new StingException("Read and Haplotype should have same length to be compared!");

        for (int i = 0; i < sz; i++) {
            Byte thisBase = this.getBase(i);
            Byte hapBase = hap.getBase(i);
            if (thisBase != null && hapBase != null && !thisBase.equals(hapBase)) {
                return false;
            }
        }
        return true;
    }
}

class PhasingStats {
    private int numReads;
    private int numVarSites;

    // Map of: sample -> PhaseCounts:
    private TreeMap<String, PhaseCounts> samplePhaseStats;

    public PhasingStats() {
        this(new TreeMap<String, PhaseCounts>());
    }

    public PhasingStats(int numReads, int numVarSites) {
        this.numReads = numReads;
        this.numVarSites = numVarSites;
        this.samplePhaseStats = new TreeMap<String, PhaseCounts>();
    }

    public PhasingStats(TreeMap<String, PhaseCounts> samplePhaseStats) {
        this.numReads = 0;
        this.numVarSites = 0;
        this.samplePhaseStats = samplePhaseStats;
    }

    public void addIn(PhasingStats other) {
        this.numReads += other.numReads;
        this.numVarSites += other.numVarSites;

        for (Map.Entry<String, PhaseCounts> sampPhaseEntry : other.samplePhaseStats.entrySet()) {
            String sample = sampPhaseEntry.getKey();
            PhaseCounts otherCounts = sampPhaseEntry.getValue();
            PhaseCounts thisCounts = this.samplePhaseStats.get(sample);
            if (thisCounts == null) {
                thisCounts = new PhaseCounts();
                this.samplePhaseStats.put(sample, thisCounts);
            }
            thisCounts.addIn(otherCounts);
        }
    }

    public int getNumReads() {
        return numReads;
    }

    public int getNumVarSites() {
        return numVarSites;
    }

    public Collection<Map.Entry<String, PhaseCounts>> getPhaseCounts() {
        return samplePhaseStats.entrySet();
    }
}

class PhaseCounts {
    public int numTestedSites; // number of het sites directly succeeding het sites
    public int numPhased;

    public PhaseCounts() {
        this.numTestedSites = 0;
        this.numPhased = 0;
    }

    public void addIn(PhaseCounts other) {
        this.numTestedSites += other.numTestedSites;
        this.numPhased += other.numPhased;
    }
}

class PhasingStatsAndOutput {
    public PhasingStats ps;
    public List<VariantContext> output;

    public PhasingStatsAndOutput(PhasingStats ps, List<VariantContext> output) {
        this.ps = ps;
        this.output = output;
    }
}