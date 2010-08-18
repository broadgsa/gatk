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

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.PreciseNonNegativeDouble;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.StingException;
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

@ReadFilters( {ZeroMappingQualityReadFilter.class} ) // Filter out all reads with zero mapping quality

public class ReadBackedPhasingWalker extends LocusWalker<PhasingStatsAndOutput, PhasingStats> {

    @Argument(fullName = "cacheWindowSize", shortName = "cacheWindow", doc = "The window size (in bases) to cache variant sites and their reads; [default:20000]", required = false)
    protected Integer cacheWindow = 20000;

    @Argument(fullName = "maxPhaseSites", shortName = "maxSites", doc = "The maximum number of successive heterozygous sites permitted to be used by the phasing algorithm; [default:20]", required = false)
    protected Integer maxPhaseSites = 20; // 2^20 == 10^6 biallelic haplotypes

    @Argument(fullName = "phaseScoreThresh", shortName = "phaseThresh", doc = "The minimum phasing quality score required to output phasing; [default:0.66]", required = false)
    protected Double phaseScoreThresh = 0.66;

    @Argument(fullName = "phasedVCFFile", shortName = "phasedVCF", doc = "The name of the phased VCF file output", required = true)
    protected String phasedVCFFile = null;

    private VCFWriter writer = null;

    private LinkedList<VariantAndReads> siteQueue = null;
    private VariantAndReads prevVr = null; // the VC emitted after phasing, and the alignment bases at the position emitted

    private static PreciseNonNegativeDouble ZERO = new PreciseNonNegativeDouble(0.0);

    private static boolean DEBUG_DETAILED = true;

    private void initializeVcfWriter(VariantContext vc) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.add(new VCFFormatHeaderLine("PQ", 1, VCFHeaderLineType.Integer, "Read-backed phasing quality"));

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
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            boolean processVariant = true;
            if (!vc.isSNP() || !vc.isBiallelic() || vc.isFiltered())
                processVariant = false;

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
                for (VariantAndReads phaseInfoVr : windowVaList)
                    logger.debug("Using phaseInfoVc = " + VariantContextUtils.getLocation(phaseInfoVr.variant));
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

                List<VariantAndReads> sampleWindowVaList = new LinkedList<VariantAndReads>();
                for (VariantAndReads phaseInfoVr : windowVaList) {
                    VariantContext phaseInfoVc = phaseInfoVr.variant;
                    Genotype phaseInfoGt = phaseInfoVc.getGenotype(samp);
                    if (phaseInfoGt.isHet()) { // otherwise, of no value to phasing
                        sampleWindowVaList.add(phaseInfoVr);
                        logger.debug("STARTING TO PHASE USING POS = " + VariantContextUtils.getLocation(phaseInfoVc));
                    }
                }
                if (sampleWindowVaList.size() > maxPhaseSites) {
                    logger.warn("Trying to phase sample " + samp + " at locus " + VariantContextUtils.getLocation(vc) + " within a window of " + cacheWindow + " bases yields " + sampleWindowVaList.size() + " heterozygous sites to phase -- REDUCING to first " + maxPhaseSites + " sites!");
                    sampleWindowVaList = sampleWindowVaList.subList(0, maxPhaseSites);
                }

                /* Will map a phase and its "complement" to a single representative phase,
                   and marginalizeTable() marginalizes to the first 2 positions [i.e., the previous position and the current position]:
                 */
                HaplotypeTableCreator tabCreator = new BiallelicComplementHaplotypeTableCreator(sampleWindowVaList, samp, 2);
                PhasingTable sampleHaps = tabCreator.getNewTable();

                // Assemble the "sub-reads" from the heterozygous positions for this sample:
                LinkedList<ReadBasesAtPosition> allPositions = new LinkedList<ReadBasesAtPosition>();
                for (VariantAndReads phaseInfoVr : sampleWindowVaList) {
                    ReadBasesAtPosition readBases = phaseInfoVr.sampleReadBases.get(samp);
                    allPositions.add(readBases);
                }
                HashMap<String, Read> allReads = convertReadBasesAtPositionToReads(allPositions);
                logger.debug("Number of reads at sites: " + allReads.size());
                int numUsedReads = 0;

                // Update the phasing table based on each of the sub-reads for this sample:
                for (Map.Entry<String, Read> nameToReads : allReads.entrySet()) {
                    Read rd = nameToReads.getValue();
                    if (rd.numNonNulls() <= 1) // can't possibly provide any phasing information, so save time
                        continue;

                    numUsedReads++;
                    if (DEBUG_DETAILED)
                        logger.debug("rd = " + rd + "\tname = " + nameToReads.getKey() + (rd.isGapped() ? "\tGAPPED" : ""));

                    for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
                        PhasingScore score = rd.matchHaplotypeClassScore(pte.getHaplotypeClass());
                        pte.getScore().integrateReadScore(score);

                        if (DEBUG_DETAILED)
                            logger.debug("score(" + rd + ", " + pte.getHaplotypeClass() + ") = " + score);
                    }
                }
                logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");
                logger.debug("numUsedReads = " + numUsedReads);

                // Marginalize each haplotype to its first 2 positions:
                sampleHaps = HaplotypeTableCreator.marginalizeTable(sampleHaps);
                logger.debug("\nPhasing table [AFTER MAPPING]:\n" + sampleHaps + "\n");

                // Determine the phase at this position:
                sampleHaps.normalizeScores();
                logger.debug("\nPhasing table [AFTER NORMALIZATION]:\n" + sampleHaps + "\n");

                PhasingTable.PhasingTableEntry maxEntry = sampleHaps.maxEntry();
                double score = maxEntry.getScore().getValue();
                logger.debug("MAX hap:\t" + maxEntry.getHaplotypeClass() + "\tscore:\t" + score);

                genotypesArePhased = (score >= phaseScoreThresh);
                if (genotypesArePhased) {
                    Biallele prevBiall = new Biallele(prevVc.getGenotype(samp));
                    ensurePhasing(biall, prevBiall, maxEntry.getHaplotypeClass().getRepresentative());
                    gtAttribs.put("PQ", new Integer(QualityUtils.probToQual(score)));

                    logger.debug("CHOSE PHASE:\n" + biall + "\n\n");
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

    public static void ensurePhasing(Biallele curBiall, Biallele prevBiall, Haplotype hap) {
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
            refBase = getSingleBase(varAllele);
        }
        else {
            refBase = vc.getReferenceBaseForIndel();
        }

        writer.add(vc, refBase);
    }

    public PhasingStats reduce(PhasingStatsAndOutput statsAndList, PhasingStats stats) {
        if (statsAndList != null) {
            writeVarContList(statsAndList.output);
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
        writeVarContList(finalList);
        if (writer != null)
            writer.close();

        out.println("Number of reads observed: " + result.getNumReads());
        out.println("Number of variant sites observed: " + result.getNumVarSites());
        out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));

        out.println("\n-- Phasing summary [minimal haplotype probability: " + phaseScoreThresh + "] --");
        for (Map.Entry<String, PhaseCounts> sampPhaseCountEntry : result.getPhaseCounts()) {
            PhaseCounts pc = sampPhaseCountEntry.getValue();
            out.println("Sample: " + sampPhaseCountEntry.getKey() + "\tNumber of tested sites: " + pc.numTestedSites + "\tNumber of phased sites: " + pc.numPhased);
        }
        out.println("");
    }

    protected void writeVarContList(List<VariantContext> varContList) {
        for (VariantContext vc : varContList) {
            writeVCF(vc);
        }
    }

    protected static HashMap<String, Read> convertReadBasesAtPositionToReads(Collection<ReadBasesAtPosition> basesAtPositions) {
        HashMap<String, Read> reads = new HashMap<String, Read>();

        int index = 0;
        for (ReadBasesAtPosition rbp : basesAtPositions) {
            for (ReadBase rb : rbp) {
                String readName = rb.readName;

                Read rd = reads.get(readName);
                if (rd == null) {
                    rd = new Read(basesAtPositions.size(), rb.mappingQual);
                    reads.put(readName, rd);
                }
                rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
            }
            index++;
        }
        return reads;
    }

    protected static byte getSingleBase(byte[] bases) {
        return bases[0];
    }

    protected static byte getSingleBase(Allele all) {
        return getSingleBase(all.getBases());
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

            return getSingleBase(topBases);
        }

        public byte getBottomBase() {
            byte[] bottomBases = bottom.getBases();
            if (bottomBases.length != 1)
                throw new StingException("LOGICAL ERROR: should not process non-SNP sites!");

            return getSingleBase(bottomBases);
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("Top:\t" + top.getBaseString() + "\n");
            sb.append("Bot:\t" + bottom.getBaseString() + "\n");
            return sb.toString();
        }

        public boolean matchesTopBase(byte base) {
            boolean matchesTop;
            if (BaseUtils.basesAreEqual(base, getTopBase()))
                matchesTop = true;
            else if (BaseUtils.basesAreEqual(base, getBottomBase()))
                matchesTop = false;
            else
                throw new StingException("LOGICAL ERROR: base MUST match either TOP or BOTTOM!");

            return matchesTop;
        }

        public byte getOtherBase(byte base) {
            byte topBase = getTopBase();
            byte botBase = getBottomBase();

            if (BaseUtils.basesAreEqual(base, topBase))
                return botBase;
            else if (BaseUtils.basesAreEqual(base, botBase))
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
                            if (!p.isDeletion()) // IGNORE deletions for now
                                readBases.putReadBase(p);
                        }
                        sampleReadBases.put(samp, readBases);
                    }
                }
            }
        }
    }

    private static class ReadBase {
        public String readName;
        public byte base;
        public int mappingQual;
        public byte baseQual;

        public ReadBase(String readName, byte base, int mappingQual, byte baseQual) {
            this.readName = readName;
            this.base = base;
            this.mappingQual = mappingQual;
            this.baseQual = baseQual;
        }
    }

    private static class ReadBasesAtPosition implements Iterable<ReadBase> {
        // list of: <read name, base>
        private LinkedList<ReadBase> bases;

        public ReadBasesAtPosition() {
            this.bases = new LinkedList<ReadBase>();
        }

        public void putReadBase(PileupElement pue) {
            ReadBase rb = new ReadBase(pue.getRead().getReadName(), pue.getBase(), pue.getMappingQual(), pue.getQual());
            bases.add(rb);
        }

        public Iterator<ReadBase> iterator() {
            return bases.iterator();
        }
    }

//
// THIS IMPLEMENTATION WILL FAIL WHEN NOT DEALING WITH SNP Alleles (e.g., MNP or INDEL), SINCE THEN THE Allele.getBases()
// FUNCTION WILL RETURN VARIABLE-LENGTH Byte ARRAYS.  IN THAT CASE, BaseArray/Haplotype/Read WILL NEED TO BE REPLACED WITH
// AN ArrayList OF Allele [OR SIMILAR OBJECT], and WON'T USE: getSingleBase(alleleI)
//
    private static abstract class HaplotypeTableCreator {
        protected Genotype[] genotypes;
        protected Biallele[] bialleles;

        public HaplotypeTableCreator(List<VariantAndReads> vaList, String sample) {
            this.genotypes = new Genotype[vaList.size()];
            this.bialleles = new Biallele[vaList.size()];

            int index = 0;
            for (VariantAndReads phaseInfoVr : vaList) {
                VariantContext phaseInfoVc = phaseInfoVr.variant;
                Genotype phaseInfoGt = phaseInfoVc.getGenotype(sample);
                genotypes[index] = phaseInfoGt;
                bialleles[index] = new Biallele(phaseInfoGt);
                index++;
            }
        }

        abstract public PhasingTable getNewTable();

        protected List<Haplotype> getAllHaplotypes() {
            int numSites = genotypes.length;
            int[] genotypeCards = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                genotypeCards[i] = genotypes[i].getPloidy();
            }

            LinkedList<Haplotype> allHaps = new LinkedList<Haplotype>();
            CardinalityCounter alleleCounter = new CardinalityCounter(genotypeCards);
            for (int[] alleleInds : alleleCounter) {
                byte[] hapBases = new byte[numSites];
                for (int i = 0; i < numSites; i++) {
                    Allele alleleI = genotypes[i].getAllele(alleleInds[i]);
                    hapBases[i] = getSingleBase(alleleI);
                }
                allHaps.add(new Haplotype(hapBases));
            }
            return allHaps;
        }

        public static PhasingTable marginalizeTable(PhasingTable table) {
            TreeMap<Haplotype, PreciseNonNegativeDouble> hapMap = new TreeMap<Haplotype, PreciseNonNegativeDouble>();
            for (PhasingTable.PhasingTableEntry pte : table) {
                Haplotype rep = pte.getHaplotypeClass().getRepresentative();
                PreciseNonNegativeDouble score = hapMap.get(rep);
                if (score == null) {
                    score = new PreciseNonNegativeDouble(0.0);
                    hapMap.put(rep, score);
                }
                score.plusEqual(pte.getScore());
            }

            PhasingTable margTable = new PhasingTable();
            for (Map.Entry<Haplotype, PreciseNonNegativeDouble> hapClassAndScore : hapMap.entrySet()) {
                Haplotype rep = hapClassAndScore.getKey();
                ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                hapList.add(rep);

                HaplotypeClass hc = new HaplotypeClass(hapList, rep);
                margTable.addEntry(hc, hapClassAndScore.getValue());
            }
            return margTable;
        }
    }

    private static class BiallelicComplementHaplotypeTableCreator extends HaplotypeTableCreator {
        private int marginalizeLength;

        public BiallelicComplementHaplotypeTableCreator(List<VariantAndReads> vaList, String sample, int marginalizeLength) {
            super(vaList, sample);
            this.marginalizeLength = marginalizeLength;
        }

        public PhasingTable getNewTable() {
            double hapClassPrior = 1.0; // can change later

            PhasingTable table = new PhasingTable();
            for (Haplotype hap : getAllHaplotypes()) {
                if (bialleles[0].matchesTopBase(hap.getBase(0))) {
                    /* hap is the "representative" haplotype [arbitrarily defined to be
                      the one with the top base at the 0th position]
                    */
                    ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                    hapList.add(hap);
                    hapList.add(complement(hap));

                    Haplotype rep = hap.subHaplotype(0, Math.min(marginalizeLength, hap.size())); // only want first marginalizeLength positions
                    HaplotypeClass hapClass = new HaplotypeClass(hapList, rep);
                    table.addEntry(hapClass, hapClassPrior);
                }
            }
            return table;
        }

        private Haplotype complement(Haplotype hap) {
            int numSites = bialleles.length;
            if (hap.size() != numSites)
                throw new StingException("INTERNAL ERROR: hap.size() != numSites");

            // Take the other base at EACH position of the Haplotype:
            byte[] complementBases = new byte[numSites];
            for (int i = 0; i < numSites; i++)
                complementBases[i] = bialleles[i].getOtherBase(hap.getBase(i));

            return new Haplotype(complementBases);
        }
    }

    private static class PhasingTable implements Iterable<PhasingTable.PhasingTableEntry> {
        private LinkedList<PhasingTableEntry> table;

        public PhasingTable() {
            this.table = new LinkedList<PhasingTableEntry>();
        }

        public PhasingTableEntry addEntry(HaplotypeClass haplotypeClass, PreciseNonNegativeDouble initialScore) {
            PhasingTableEntry pte = new PhasingTableEntry(haplotypeClass, new PhasingScore(initialScore));
            table.add(pte);
            return pte;
        }

        public PhasingTableEntry addEntry(HaplotypeClass haplotypeClass, double initialScore) {
            return addEntry(haplotypeClass, new PreciseNonNegativeDouble(initialScore));
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
                if (maxPte == null || pte.getScore().gt(maxPte.getScore())) {
                    maxPte = pte;
                }
            }
            return maxPte;
        }

        public void normalizeScores() {
            PreciseNonNegativeDouble normalizeBy = new PreciseNonNegativeDouble(ZERO);
            for (PhasingTableEntry pte : table)
                normalizeBy.plusEqual(pte.getScore());

            if (!normalizeBy.equals(ZERO)) { // prevent precision problems
                logger.debug("normalizeBy = " + normalizeBy);
                for (PhasingTableEntry pte : table)
                    pte.getScore().divEqual(normalizeBy);
            }
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("-------------------\n");
            for (PhasingTableEntry pte : this) {
                sb.append("Haplotypes:\t" + pte.getHaplotypeClass() + "\tScore:\t" + pte.getScore() + "\n");
            }
            sb.append("-------------------\n");
            return sb.toString();
        }

        public static class PhasingTableEntry implements Comparable<PhasingTableEntry> {
            private HaplotypeClass haplotypeClass;
            private PhasingScore score;

            public PhasingTableEntry(HaplotypeClass haplotypeClass, PhasingScore score) {
                this.haplotypeClass = haplotypeClass;
                this.score = score;
            }

            public HaplotypeClass getHaplotypeClass() {
                return haplotypeClass;
            }

            public PhasingScore getScore() {
                return score;
            }

            public int compareTo(PhasingTableEntry that) {
                return this.getScore().compareTo(that.getScore());
            }
        }
    }
}



class PhasingScore extends PreciseNonNegativeDouble {
    public PhasingScore(double score) {
        super(score);
    }

    public PhasingScore(PreciseNonNegativeDouble val) {
        super(val);
    }

    public PhasingScore integrateReadScore(PhasingScore score) {
        timesEqual(score);
        return this;
    }
}

abstract class BaseArray implements Comparable<BaseArray> {
    protected Byte[] bases;

    public BaseArray(byte[] bases) {
        this.bases = new Byte[bases.length];
        for (int i = 0; i < bases.length; i++)
            this.bases[i] = bases[i];
    }

    public BaseArray(Byte[] bases) {
        this.bases = Arrays.copyOf(bases, bases.length);
    }

    public BaseArray(int length) {
        this.bases = new Byte[length];
        Arrays.fill(this.bases, null);
    }

    public BaseArray(BaseArray other) {
        this(other.bases);
    }

    public void updateBase(int index, Byte base) {
        bases[index] = base;
    }

    public Byte getBase(int index) {
        return bases[index];
    }

    public int size() {
        return bases.length;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(bases.length);
        for (Byte b : bases)
            sb.append(b != null ? (char) b.byteValue() : "_");

        return sb.toString();
    }

    public int compareTo(BaseArray that) {
        int sz = this.bases.length;
        if (sz != that.bases.length)
            return (sz - that.bases.length);

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
            else if (!BaseUtils.basesAreEqual(thisBase, thatBase)) {
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

    public Haplotype(Haplotype other) {
        super(other);
    }

    public void updateBase(int index, Byte base) {
        if (base == null) {
            throw new StingException("Internal error: CANNOT have null for a missing Haplotype base!");
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
        return new Haplotype(Arrays.copyOfRange(bases, fromIndex, toIndex));
    }
}

class HaplotypeClass implements Iterable<Haplotype> {
    private ArrayList<Haplotype> haps;
    private Haplotype rep;

    public HaplotypeClass(ArrayList<Haplotype> haps, Haplotype rep) {
        this.haps = haps;
        this.rep = rep;
    }

    public Iterator<Haplotype> iterator() {
        return haps.iterator();
    }

    public Haplotype getRepresentative() {
        return rep;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        boolean isFirst = true;
        for (Haplotype h : haps) {
            if (isFirst)
                isFirst = false;
            else
                sb.append("|");

            sb.append(h);
        }
        sb.append(" [").append(rep).append("]");
        return sb.toString();
    }
}

class Read extends BaseArray {
    private PreciseNonNegativeDouble mappingProb; // the probability that this read is mapped correctly
    private PreciseNonNegativeDouble[] baseProbs; // the probabilities that the base identities are CORRECT
    private PreciseNonNegativeDouble[] baseErrorProbs; // the probabilities that the base identities are INCORRECT

    public Read(int length, int mappingQual) {
        super(length);

        this.mappingProb = new PreciseNonNegativeDouble(QualityUtils.qualToProb(mappingQual));

        this.baseProbs = new PreciseNonNegativeDouble[length];
        Arrays.fill(this.baseProbs, null);

        this.baseErrorProbs = new PreciseNonNegativeDouble[length];
        Arrays.fill(this.baseErrorProbs, null);
    }

    public void updateBaseAndQuality(int index, Byte base, byte baseQual) {
        updateBase(index, base);

        double errProb = QualityUtils.qualToErrorProb(baseQual);
        baseProbs[index] = new PreciseNonNegativeDouble(1.0 - errProb);
        baseErrorProbs[index] = new PreciseNonNegativeDouble(errProb);
    }

    public int numNonNulls() {
        int num = 0;
        for (int i = 0; i < bases.length; i++) {
            if (getBase(i) != null)
                num++;
        }
        return num;
    }

    public PhasingScore matchHaplotypeClassScore(HaplotypeClass hapClass) {
        PreciseNonNegativeDouble value = new PreciseNonNegativeDouble(0.0);
        for (Haplotype h : hapClass)
            value.plusEqual(matchHaplotypeScore(h));

        return new PhasingScore(value);
    }

    private PreciseNonNegativeDouble matchHaplotypeScore(Haplotype hap) {
        PreciseNonNegativeDouble score = new PreciseNonNegativeDouble(1.0);

        int sz = this.bases.length;
        if (sz != hap.bases.length)
            throw new StingException("Read and Haplotype should have same length to be compared!");

        for (int i = 0; i < sz; i++) {
            Byte thisBase = this.getBase(i);
            Byte hapBase = hap.getBase(i);
            if (thisBase != null && hapBase != null) {
                if (BaseUtils.basesAreEqual(thisBase, hapBase))
                    score.timesEqual(baseProbs[i]);
                else
                    score.timesEqual(baseErrorProbs[i]);
                score.timesEqual(mappingProb);
            }
        }
        return score;
    }

    private enum ReadStage {
        BEFORE_BASES, BASES_1, NO_BASES, BASES_2
    }

    public boolean isGapped() {
        ReadStage s = ReadStage.BEFORE_BASES;

        for (int i = 0; i < bases.length; i++) {
            if (getBase(i) != null) { // has a base at i
                if (s == ReadStage.BEFORE_BASES)
                    s = ReadStage.BASES_1;
                else if (s == ReadStage.NO_BASES) {
                    s = ReadStage.BASES_2;
                    break;
                }
            }
            else { // no base at i
                if (s == ReadStage.BASES_1)
                    s = ReadStage.NO_BASES;
            }
        }

        return (s == ReadStage.BASES_2);
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

class CardinalityCounter implements Iterator<int[]>, Iterable<int[]> {
    private int[] cards;
    private int[] valList;
    private boolean hasNext;

    public CardinalityCounter(int[] cards) {
        this.cards = cards;
        this.valList = new int[cards.length];
        for (int i = 0; i < cards.length; i++) {
            if (this.cards[i] <= 0)
                throw new StingException("CANNOT have zero cardinalities!");
            this.valList[i] = 0;
        }
        this.hasNext = true;
    }

    public boolean hasNext() {
        return hasNext;
    }

    public int[] next() {
        if (!hasNext())
            throw new StingException("CANNOT iterate past end!");

        // Copy the assignment to be returned:
        int[] nextList = new int[valList.length];
        for (int i = 0; i < valList.length; i++)
            nextList[i] = valList[i];

        // Find the assignment after this one:
        hasNext = false;
        int i = cards.length - 1;
        for (; i >= 0; i--) {
            if (valList[i] < (cards[i] - 1)) {
                valList[i]++;
                hasNext = true;
                break;
            }
            valList[i] = 0;
        }

        return nextList;
    }

    public void remove() {
        throw new StingException("Cannot remove from CardinalityCounter!");
    }

    public Iterator<int[]> iterator() {
        return this;
    }

}
