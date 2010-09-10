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
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.playground.gatk.walkers.phasing.*;

import java.io.*;
import java.util.*;

import static org.broadinstitute.sting.utils.vcf.VCFUtils.getVCFHeadersFromRods;


/**
 * Walks along all variant ROD loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using upstream and downstream reads).
 */
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE}, referenceMetaData = @RMD(name = "variant", type = ReferenceOrderedDatum.class))
@By(DataSource.READS)

@ReadFilters({ZeroMappingQualityReadFilter.class})
// Filter out all reads with zero mapping quality

public class ReadBackedPhasingWalker extends RodWalker<PhasingStatsAndOutput, PhasingStats> {

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;

    @Argument(fullName = "cacheWindowSize", shortName = "cacheWindow", doc = "The window size (in bases) to cache variant sites and their reads; [default:20000]", required = false)
    protected Integer cacheWindow = 20000;

    @Argument(fullName = "maxPhaseSites", shortName = "maxSites", doc = "The maximum number of successive heterozygous sites permitted to be used by the phasing algorithm; [default:10]", required = false)
    protected Integer maxPhaseSites = 10; // 2^10 == 10^3 biallelic haplotypes

    @Argument(fullName = "phaseQualityThresh", shortName = "phaseThresh", doc = "The minimum phasing quality score required to output phasing; [default:10.0]", required = false)
    protected Double phaseQualityThresh = 10.0; // PQ = 10.0 <=> P(error) = 10^(-10/10) = 0.1, P(correct) = 0.9

    @Argument(fullName = "variantStatsFilePrefix", shortName = "variantStats", doc = "The prefix of the VCF/phasing statistics files", required = false)
    protected String variantStatsFilePrefix = null;

    private LinkedList<VariantAndReads> unphasedSiteQueue = null;
    private LinkedList<VariantAndReads> phasedSites = null; // the phased VCs to be emitted, and the alignment bases at these positions

    private static PreciseNonNegativeDouble ZERO = new PreciseNonNegativeDouble(0.0);

    private LinkedList<String> rodNames = null;
    private PhasingQualityStatsWriter statsWriter = null;

    public void initialize() {
        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        unphasedSiteQueue = new LinkedList<VariantAndReads>();
        phasedSites = new LinkedList<VariantAndReads>();

        initializeVcfWriter();

        if (variantStatsFilePrefix != null)
            statsWriter = new PhasingQualityStatsWriter(variantStatsFilePrefix);
    }

    private void initializeVcfWriter() {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.add(new VCFFormatHeaderLine("PQ", 1, VCFHeaderLineType.Float, "Read-backed phasing quality"));

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), rodNames);
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(rodNameToHeader.get(rodNames.get(0)).getGenotypeSamples())));
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public PhasingStats reduceInit() {
        return new PhasingStats();
    }

    /**
     * For each site of interest, cache the current site and then use the cache to phase all sites
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

        boolean requireStartHere = true; // only see each VariantContext once
        boolean takeFirstOnly = false; // take as many entries as the VCF file has
        for (VariantContext vc : tracker.getVariantContexts(ref, rodNames, null, context.getLocation(), requireStartHere, takeFirstOnly)) {
            boolean processVariant = true;
            if (!isUnfilteredBiallelicSNP(vc))
                processVariant = false;

            VariantAndReads vr = new VariantAndReads(vc, context, processVariant);
            unphasedSiteQueue.add(vr);
            logger.debug("Added variant to queue = " + VariantContextUtils.getLocation(vr.variant));

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
        List<VariantContext> phasedList = processQueue(phaseStats, false);

        return new PhasingStatsAndOutput(phaseStats, phasedList);
    }

    private List<VariantContext> processQueue(PhasingStats phaseStats, boolean processAll) {
        List<VariantContext> oldPhasedList = new LinkedList<VariantContext>();

        if (!unphasedSiteQueue.isEmpty()) {
            GenomeLoc lastLocus = null;
            if (!processAll)
                lastLocus = VariantContextUtils.getLocation(unphasedSiteQueue.peekLast().variant);

            while (!unphasedSiteQueue.isEmpty()) {
                if (!processAll) { // otherwise, phase until the end of unphasedSiteQueue
                    VariantContext nextToPhaseVc = unphasedSiteQueue.peek().variant;
                    if (isInWindowRange(lastLocus, VariantContextUtils.getLocation(nextToPhaseVc))) {
                        /* lastLocus is still not far enough ahead of nextToPhaseVc to have all phasing information for nextToPhaseVc
                          (note that we ASSUME that the VCF is ordered by <contig,locus>) */
                        break;
                    }
                    // Already saw all variant positions within cacheWindow distance ahead of vc (on its contig)
                }
                // Update phasedSites before it's used in finalizePhasing:
                oldPhasedList.addAll(discardIrrelevantPhasedSites());
                logger.debug("oldPhasedList(1) = " + toStringVCL(oldPhasedList));

                VariantAndReads phasedVr = finalizePhasing(unphasedSiteQueue.remove(), phaseStats);
                logger.debug("Finalized phasing for " + VariantContextUtils.getLocation(phasedVr.variant));
                phasedSites.add(phasedVr);
            }
        }

        // Update phasedSites after finalizePhasing is done:
        oldPhasedList.addAll(discardIrrelevantPhasedSites());
        logger.debug("oldPhasedList(2) = " + toStringVCL(oldPhasedList));
        return oldPhasedList;
    }

    private List<VariantContext> discardIrrelevantPhasedSites() {
        List<VariantContext> vcList = new LinkedList<VariantContext>();

        VariantContext nextToPhaseVc = null;
        if (!unphasedSiteQueue.isEmpty())
            nextToPhaseVc = unphasedSiteQueue.peek().variant;

        while (!phasedSites.isEmpty()) {
            VariantAndReads phasedVr = phasedSites.peek();
            VariantContext phasedVc = phasedVr.variant;
            if (nextToPhaseVc != null && phasedVr.processVariant && isInWindowRange(phasedVc, nextToPhaseVc)) {
                // nextToPhaseVc is still not far enough ahead of phasedVc to exclude phasedVc from calculations
                break;
            }
            vcList.add(phasedSites.remove().variant);
        }

        return vcList;
    }

    /* Phase vc (removed head of unphasedSiteQueue) using all VariantContext objects in
       phasedSites, and all in unphasedSiteQueue that are within cacheWindow distance ahead of vc (on its contig).

       ASSUMES: All VariantContexts in unphasedSiteQueue are in positions downstream of vc (head of queue).
     */

    private VariantAndReads finalizePhasing(VariantAndReads vr, PhasingStats phaseStats) {
        if (!vr.processVariant)
            return vr; // return vr as is

        // Find the previous VariantContext (that was processed and phased):
        VariantAndReads prevVr = null;
        Iterator<VariantAndReads> backwardsIt = phasedSites.descendingIterator(); // look at most recently phased sites
        while (backwardsIt.hasNext()) {
            VariantAndReads backVr = backwardsIt.next();
            if (backVr.processVariant) {
                prevVr = backVr;
                break;
            }
        }
        if (prevVr == null)
            return vr; // return vr as is, since cannot phase against "nothing" (vc is at the beginning of the chromosome, or the previous was so far back it was removed from phasedSites)

        VariantContext vc = vr.variant;
        logger.debug("Will phase vc = " + VariantContextUtils.getLocation(vc));

        LinkedList<VariantAndReads> windowVaList = new LinkedList<VariantAndReads>();

        // Include previously phased sites in the phasing computation:
        for (VariantAndReads phasedVr : phasedSites) {
            if (phasedVr.processVariant)
                windowVaList.add(phasedVr);
        }

        // Add position to be phased:
        windowVaList.add(vr);

        // Include as of yet unphased sites in the phasing computation:
        for (VariantAndReads nextVr : unphasedSiteQueue) {
            if (!isInWindowRange(vc, nextVr.variant)) //nextVr too far ahead of the range used for phasing vc
                break;
            if (nextVr.processVariant) // include in the phasing computation
                windowVaList.add(nextVr);
        }

        if (logger.isDebugEnabled()) {
            for (VariantAndReads phaseInfoVr : windowVaList)
                logger.debug("Using phaseInfoVc = " + VariantContextUtils.getLocation(phaseInfoVr.variant));
        }
        logger.debug("");

        Map<String, Genotype> sampGenotypes = vc.getGenotypes();
        Map<String, Genotype> phasedGtMap = new TreeMap<String, Genotype>();

        // Perform per-sample phasing:
        TreeMap<String, PhaseCounts> samplePhaseStats = new TreeMap<String, PhaseCounts>();
        for (Map.Entry<String, Genotype> sampGtEntry : sampGenotypes.entrySet()) {
            logger.debug("sample = " + sampGtEntry.getKey());
            boolean genotypesArePhased = true; // phase by default

            String samp = sampGtEntry.getKey();
            Genotype gt = sampGtEntry.getValue();
            BialleleSNP biall = new BialleleSNP(gt);
            HashMap<String, Object> gtAttribs = new HashMap<String, Object>(gt.getAttributes());

            if (gt.isHet()) {
                VariantContext prevVc = prevVr.variant;
                Genotype prevGenotype = prevVc.getGenotype(samp);
                if (prevGenotype.isHet()) { //otherwise, can trivially phase
                    logger.debug("NON-TRIVIALLY CARE about TOP vs. BOTTOM for: " + "\n" + biall);

                    List<VariantAndReads> sampleWindowVaList = new LinkedList<VariantAndReads>();
                    int phasingSiteIndex = -1;
                    int currentIndex = 0;
                    for (VariantAndReads phaseInfoVr : windowVaList) {
                        VariantContext phaseInfoVc = phaseInfoVr.variant;
                        Genotype phaseInfoGt = phaseInfoVc.getGenotype(samp);
                        if (phaseInfoGt.isHet()) { // otherwise, of no value to phasing
                            sampleWindowVaList.add(phaseInfoVr);
                            if (phasingSiteIndex == -1) {
                                if (phaseInfoVr == vr)
                                    phasingSiteIndex = currentIndex; // index of vr in sampleWindowVaList
                                else
                                    currentIndex++;
                            }
                            logger.debug("STARTING TO PHASE USING POS = " + VariantContextUtils.getLocation(phaseInfoVc));
                        }
                    }
                    if (logger.isDebugEnabled() && (phasingSiteIndex == -1 || phasingSiteIndex == 0))
                        throw new GATKException("Internal error: could NOT find vr and/or prevVr!");

                    if (sampleWindowVaList.size() > maxPhaseSites) {
                        logger.warn("Trying to phase sample " + samp + " at locus " + VariantContextUtils.getLocation(vc) + " within a window of " + cacheWindow + " bases yields " + sampleWindowVaList.size() + " heterozygous sites to phase:\n" + toStringVRL(sampleWindowVaList));

                        int prevSiteIndex = phasingSiteIndex - 1; // index of prevVr in sampleWindowVaList
                        int numToUse = maxPhaseSites - 2; // since always keep prevVr and vr

                        int numOnLeft = prevSiteIndex;
                        int numOnRight = sampleWindowVaList.size() - (phasingSiteIndex + 1);

                        int useOnLeft, useOnRight;
                        if (numOnLeft <= numOnRight) {
                            int halfToUse = new Double(Math.floor(numToUse / 2.0)).intValue(); // skimp on the left [floor], and be generous with the right side
                            useOnLeft = Math.min(halfToUse, numOnLeft);
                            useOnRight = Math.min(numToUse - useOnLeft, numOnRight);
                        }
                        else { // numOnRight < numOnLeft
                            int halfToUse = new Double(Math.ceil(numToUse / 2.0)).intValue(); // be generous with the right side [ceil]
                            useOnRight = Math.min(halfToUse, numOnRight);
                            useOnLeft = Math.min(numToUse - useOnRight, numOnLeft);
                        }
                        int startIndex = prevSiteIndex - useOnLeft;
                        int stopIndex = phasingSiteIndex + useOnRight + 1; // put the index 1 past the desired index to keep
                        phasingSiteIndex -= startIndex;
                        sampleWindowVaList = sampleWindowVaList.subList(startIndex, stopIndex);
                        logger.warn("REDUCED to " + sampleWindowVaList.size() + " sites:\n" + toStringVRL(sampleWindowVaList));
                    }

                    PhaseResult pr = phaseSample(samp, sampleWindowVaList, phasingSiteIndex);
                    genotypesArePhased = (pr.phaseQuality >= phaseQualityThresh);
                    if (genotypesArePhased) {
                        BialleleSNP prevBiall = new BialleleSNP(prevGenotype);

                        logger.debug("THE PHASE PREVIOUSLY CHOSEN FOR PREVIOUS:\n" + prevBiall + "\n");
                        logger.debug("THE PHASE CHOSEN HERE:\n" + biall + "\n\n");

                        ensurePhasing(biall, prevBiall, pr.haplotype);
                        gtAttribs.put("PQ", pr.phaseQuality);
                    }

                    if (statsWriter != null)
                        statsWriter.addStat(samp, VariantContextUtils.getLocation(vc), distance(prevVc, vc), pr.phaseQuality, pr.numReads);

                    PhaseCounts sampPhaseCounts = samplePhaseStats.get(samp);
                    if (sampPhaseCounts == null) {
                        sampPhaseCounts = new PhaseCounts();
                        samplePhaseStats.put(samp, sampPhaseCounts);
                    }
                    sampPhaseCounts.numTestedSites++;
                    if (genotypesArePhased)
                        sampPhaseCounts.numPhased++;
                }
            }
            List<Allele> phasedAll = biall.getAllelesAsList();
            Genotype phasedGt = new Genotype(gt.getSampleName(), phasedAll, gt.getNegLog10PError(), gt.getFilters(), gtAttribs, genotypesArePhased);
            phasedGtMap.put(samp, phasedGt);
        }
        phaseStats.addIn(new PhasingStats(samplePhaseStats));

        VariantContext phasedVc = new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), phasedGtMap, vc.getNegLog10PError(), vc.getFilters(), vc.getAttributes());
        return new VariantAndReads(phasedVc, vr.sampleReadBases, vr.processVariant);
    }

    private PhaseResult phaseSample(String sample, List<VariantAndReads> variantList, int phasingSiteIndex) {
        /* Will map a phase and its "complement" to a single representative phase,
          and marginalizeTable() marginalizes to 2 positions [starting at the previous position, and then the current position]:
        */
        HaplotypeTableCreator tabCreator = new BiallelicComplementHaplotypeTableCreator(variantList, sample, phasingSiteIndex - 1, 2);
        PhasingTable sampleHaps = tabCreator.getNewTable();

        // Assemble the "sub-reads" from the heterozygous positions for this sample:
        LinkedList<ReadBasesAtPosition> allPositions = new LinkedList<ReadBasesAtPosition>();
        for (VariantAndReads phaseInfoVr : variantList) {
            ReadBasesAtPosition readBases = phaseInfoVr.sampleReadBases.get(sample);
            if (readBases == null)
                readBases = new ReadBasesAtPosition(); // for transparency, put an empty list of bases at this position for sample
            allPositions.add(readBases);
        }
        HashMap<String, Read> allReads = convertReadBasesAtPositionToReads(allPositions);
        logger.debug("Number of TOTAL reads [including those covering only 1 position] at sites: " + allReads.size());
        int numUsedReads = 0;

        // Update the phasing table based on each of the sub-reads for this sample:
        for (Map.Entry<String, Read> nameToReads : allReads.entrySet()) {
            Read rd = nameToReads.getValue();
            if (rd.numNonNulls() <= 1) // can't possibly provide any phasing information, so save time
                continue;

            numUsedReads++;
            logger.debug("rd = " + rd + "\tname = " + nameToReads.getKey() + (rd.isGapped() ? "\tGAPPED" : ""));

            for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
                PhasingScore score = rd.matchHaplotypeClassScore(pte.getHaplotypeClass());
                pte.getScore().integrateReadScore(score);

                logger.debug("score(" + rd + ", " + pte.getHaplotypeClass() + ") = " + score);
            }
        }
        logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");
        logger.debug("numUsedReads [covering > 1 position in the haplotype] = " + numUsedReads);

        // Marginalize each haplotype to its first 2 positions:
        sampleHaps = HaplotypeTableCreator.marginalizeTable(sampleHaps);
        logger.debug("\nPhasing table [AFTER MAPPING]:\n" + sampleHaps + "\n");

        // Determine the phase at this position:
        sampleHaps.normalizeScores();
        logger.debug("\nPhasing table [AFTER NORMALIZATION]:\n" + sampleHaps + "\n");

        PhasingTable.PhasingTableEntry maxEntry = sampleHaps.maxEntry();
        double posteriorProb = maxEntry.getScore().getValue();

        // convert posteriorProb to PHRED scale, but do NOT cap the quality as in QualityUtils.probToQual(posteriorProb):
        PreciseNonNegativeDouble sumErrorProbs = new PreciseNonNegativeDouble(ZERO);
        for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
            if (pte != maxEntry)
                sumErrorProbs.plusEqual(pte.getScore());
        }
        double phaseQuality = -10.0 * (sumErrorProbs.getLog10Value());

        logger.debug("MAX hap:\t" + maxEntry.getHaplotypeClass() + "\tposteriorProb:\t" + posteriorProb + "\tphaseQuality:\t" + phaseQuality);

        return new PhaseResult(maxEntry.getHaplotypeClass().getRepresentative(), phaseQuality, numUsedReads);
    }

    /*
        Ensure that curBiall is phased relative to prevBiall as specified by hap.
     */

    public static void ensurePhasing(BialleleSNP curBiall, BialleleSNP prevBiall, Haplotype hap) {
        if (hap.size() < 2)
            throw new GATKException("LOGICAL ERROR: Only considering haplotypes of length > 2!");

        byte prevBase = hap.getBase(0); // The 1st base in the haplotype
        byte curBase = hap.getBase(1);  // The 2nd base in the haplotype

        boolean chosePrevTopChrom = prevBiall.matchesTopBase(prevBase);
        boolean choseCurTopChrom = curBiall.matchesTopBase(curBase);
        if (chosePrevTopChrom != choseCurTopChrom)
            curBiall.swapAlleles();
    }

    private boolean isInWindowRange(VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(vc2);

        return isInWindowRange(loc1, loc2);
    }

    private boolean isInWindowRange(GenomeLoc loc1, GenomeLoc loc2) {
        return (loc1.onSameContig(loc2) && loc1.distance(loc2) <= cacheWindow);
    }

    private static int distance(VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(vc2);
        if (!loc1.onSameContig(loc2))
            return Integer.MAX_VALUE;

        return loc1.distance(loc2);
    }

    private void writeVCF(VariantContext vc) {
        byte refBase;
        if (!vc.isIndel()) {
            Allele varAllele = vc.getReference();
            refBase = BialleleSNP.getSingleBase(varAllele);
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
     * Phase anything left in the cached unphasedSiteQueue, and report the number of reads and VariantContexts processed.
     *
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(PhasingStats result) {
        List<VariantContext> finalList = processQueue(result, true);
        writeVarContList(finalList);
        if (statsWriter != null)
            statsWriter.close();

        System.out.println("Number of reads observed: " + result.getNumReads());
        System.out.println("Number of variant sites observed: " + result.getNumVarSites());
        System.out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));

        System.out.println("\n-- Phasing summary [minimal haplotype probability: " + phaseQualityThresh + "] --");
        for (Map.Entry<String, PhaseCounts> sampPhaseCountEntry : result.getPhaseCounts()) {
            PhaseCounts pc = sampPhaseCountEntry.getValue();
            System.out.println("Sample: " + sampPhaseCountEntry.getKey() + "\tNumber of tested sites: " + pc.numTestedSites + "\tNumber of phased sites: " + pc.numPhased);
        }
        System.out.println("");
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


    /*
       Inner classes:
     */
    private static class VariantAndReads {
        public VariantContext variant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;
        public boolean processVariant;

        public VariantAndReads(VariantContext variant, HashMap<String, ReadBasesAtPosition> sampleReadBases, boolean processVariant) {
            this.variant = variant;
            this.sampleReadBases = sampleReadBases;
            this.processVariant = processVariant;
        }

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

    private static String toStringVRL(List<VariantAndReads> vrList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (VariantAndReads vr : vrList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(VariantContextUtils.getLocation(vr.variant));
        }
        return sb.toString();
    }

    private static String toStringVCL(List<VariantContext> vcList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (VariantContext vc : vcList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(VariantContextUtils.getLocation(vc));
        }
        return sb.toString();
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

        public HaplotypeTableCreator(List<VariantAndReads> vaList, String sample) {
            this.genotypes = new Genotype[vaList.size()];
            int index = 0;
            for (VariantAndReads phaseInfoVr : vaList) {
                VariantContext phaseInfoVc = phaseInfoVr.variant;
                Genotype phaseInfoGt = phaseInfoVc.getGenotype(sample);
                genotypes[index++] = phaseInfoGt;
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
                    hapBases[i] = BialleleSNP.getSingleBase(alleleI);
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
                    score = new PreciseNonNegativeDouble(ZERO);
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
        private BialleleSNP[] bialleleSNPs;
        private int startIndex;
        private int marginalizeLength;

        public BiallelicComplementHaplotypeTableCreator(List<VariantAndReads> vaList, String sample, int startIndex, int marginalizeLength) {
            super(vaList, sample);

            this.bialleleSNPs = new BialleleSNP[genotypes.length];
            for (int i = 0; i < genotypes.length; i++)
                bialleleSNPs[i] = new BialleleSNP(genotypes[i]);

            this.startIndex = startIndex;
            this.marginalizeLength = marginalizeLength;
        }

        public PhasingTable getNewTable() {
            double hapClassPrior = 1.0; // can change later

            PhasingTable table = new PhasingTable();
            for (Haplotype hap : getAllHaplotypes()) {
                if (bialleleSNPs[startIndex].matchesTopBase(hap.getBase(startIndex))) {
                    /* hap is the "representative" haplotype [DEFINED here to be
                      the one with the top base at the startIndex position.
                      NOTE that it is CRITICAL that this definition be consistent with the representative sub-haplotypes defined below!]
                    */
                    ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                    hapList.add(hap);
                    hapList.add(complement(hap));

                    // want marginalizeLength positions starting at startIndex:
                    Haplotype rep = hap.subHaplotype(startIndex, startIndex + marginalizeLength);

                    HaplotypeClass hapClass = new HaplotypeClass(hapList, rep);
                    table.addEntry(hapClass, hapClassPrior);
                }
            }
            return table;
        }

        private Haplotype complement(Haplotype hap) {
            int numSites = bialleleSNPs.length;
            if (hap.size() != numSites)
                throw new GATKException("INTERNAL ERROR: hap.size() != numSites");

            // Take the other base at EACH position of the Haplotype:
            byte[] complementBases = new byte[numSites];
            for (int i = 0; i < numSites; i++)
                complementBases[i] = bialleleSNPs[i].getOtherBase(hap.getBase(i));

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

    private static class PhaseResult {
        public Haplotype haplotype;
        public double phaseQuality;
        public int numReads;

        public PhaseResult(Haplotype haplotype, double phaseQuality, int numReads) {
            this.haplotype = haplotype;
            this.phaseQuality = phaseQuality;
            this.numReads = numReads;
        }
    }

    public static boolean isUnfilteredBiallelicSNP(VariantContext vc) {
        return (vc.isSNP() && vc.isBiallelic() && !vc.isFiltered());
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
            throw new GATKException("Internal error: CANNOT have null for a missing Haplotype base!");
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
        return new Haplotype(Arrays.copyOfRange(bases, fromIndex, Math.min(toIndex, size())));
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
            throw new GATKException("Read and Haplotype should have same length to be compared!");

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

class PhasingQualityStatsWriter {
    private String variantStatsFilePrefix;
    private HashMap<String, BufferedWriter> sampleToStatsWriter = new HashMap<String, BufferedWriter>();

    public PhasingQualityStatsWriter(String variantStatsFilePrefix) {
        this.variantStatsFilePrefix = variantStatsFilePrefix;
    }

    public void addStat(String sample, GenomeLoc locus, int distanceFromPrevious, double phasingQuality, int numReads) {
        BufferedWriter sampWriter = sampleToStatsWriter.get(sample);
        if (sampWriter == null) {
            String fileName = variantStatsFilePrefix + "." + sample + ".locus_distance_PQ_numReads.txt";

            FileOutputStream output;
            try {
                output = new FileOutputStream(fileName);
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Unable to create phasing quality stats file at location: " + fileName);
            }
            sampWriter = new BufferedWriter(new OutputStreamWriter(output));
            sampleToStatsWriter.put(sample, sampWriter);
        }
        try {
            sampWriter.write(locus + "\t" + distanceFromPrevious + "\t" + phasingQuality + "\t" + numReads + "\n");
            sampWriter.flush();
        } catch (IOException e) {
            throw new RuntimeException("Unable to write to per-sample phasing quality stats file", e);
        }
    }

    public void close() {
        for (Map.Entry<String, BufferedWriter> sampWriterEntry : sampleToStatsWriter.entrySet()) {
            BufferedWriter sampWriter = sampWriterEntry.getValue();
            try {
                sampWriter.flush();
                sampWriter.close();
            } catch (IOException e) {
                throw new RuntimeException("Unable to close per-sample phasing quality stats file");
            }
        }
    }
}
