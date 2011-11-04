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
package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.MappingQualityZeroFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.DisjointSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.*;
import java.util.*;

import static org.broadinstitute.sting.utils.codecs.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using upstream and downstream reads).
 *
 * <p>
 * Performs physical phasing of SNP calls, based on sequencing reads.
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * VCF file of SNP calls, BAM file of sequence reads.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Phased VCF file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ReadBackedPhasing
 *      -R reference.fasta
 *      -I reads.bam
 *      --variant SNPs.vcf
 *      -L SNPs.vcf
 *      -o phased_SNPs.vcf
 *      --phaseQualityThresh 20.0
 * </pre>
 *
 * @author Menachem Fromer
 * @since July 2010
 */
@Allows(value = {DataSource.READS, DataSource.REFERENCE})
@Requires(value = {DataSource.READS, DataSource.REFERENCE})
@By(DataSource.READS)

// Filter out all reads with zero mapping quality
@ReadFilters({MappingQualityZeroFilter.class})

public class ReadBackedPhasingWalker extends RodWalker<PhasingStatsAndOutput, PhasingStats> {
    private static final boolean DEBUG = false;
    /**
     * The VCF file we are phasing variants from.
     *
     * All heterozygous variants found in this VCF file will be phased, where possible
     */
    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc = "File to which variants should be written", required = true)
    protected VCFWriter writer = null;

    @Argument(fullName = "cacheWindowSize", shortName = "cacheWindow", doc = "The window size (in bases) to cache variant sites and their reads for the phasing procedure", required = false)
    protected Integer cacheWindow = 20000;

    @Argument(fullName = "maxPhaseSites", shortName = "maxSites", doc = "The maximum number of successive heterozygous sites permitted to be used by the phasing algorithm", required = false)
    protected Integer maxPhaseSites = 10; // 2^10 == 10^3 diploid haplotypes

    @Argument(fullName = "phaseQualityThresh", shortName = "phaseThresh", doc = "The minimum phasing quality score required to output phasing", required = false)
    protected Double phaseQualityThresh = 10.0; // PQ = 10.0 <=> P(error) = 10^(-10/10) = 0.1, P(correct) = 0.9

    @Hidden
    @Argument(fullName = "variantStatsFilePrefix", shortName = "variantStats", doc = "The prefix of the VCF/phasing statistics files [For DEBUGGING purposes only - DO NOT USE!]", required = false)
    protected String variantStatsFilePrefix = null;
    private PhasingQualityStatsWriter statsWriter = null;

    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for phasing", required = false)
    public int MIN_BASE_QUALITY_SCORE = 17;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for phasing", required = false)
    public int MIN_MAPPING_QUALITY_SCORE = 20;

    @Argument(fullName = "sampleToPhase", shortName = "sampleToPhase", doc = "Only include these samples when phasing", required = false)
    protected List<String> samplesToPhase = null;

    private GenomeLoc mostDownstreamLocusReached = null;

    private LinkedList<VariantAndReads> unphasedSiteQueue = null;
    private CloneableIteratorLinkedList<UnfinishedVariantAndReads> partiallyPhasedSites = null; // the phased VCs to be emitted, and the alignment bases at these positions

    private static PreciseNonNegativeDouble ZERO = new PreciseNonNegativeDouble(0.0);

    public static final String PQ_KEY = "PQ";

    // In order to detect phase inconsistencies:
    private static final double FRACTION_OF_MEAN_PQ_CHANGES = 0.1; // If the PQ decreases by this fraction of the mean PQ changes (thus far), then this read is inconsistent with previous reads
    private static final double MAX_FRACTION_OF_INCONSISTENT_READS = 0.1; // If there are more than this fraction of inconsistent reads, then flag this site

    public static final String PHASING_INCONSISTENT_KEY = "PhasingInconsistent";

    @Argument(fullName = "enableMergePhasedSegregatingPolymorphismsToMNP", shortName = "enableMergeToMNP", doc = "Merge consecutive phased sites into MNP records", required = false)
    protected boolean enableMergePhasedSegregatingPolymorphismsToMNP = false;

    @Argument(fullName = "maxGenomicDistanceForMNP", shortName = "maxDistMNP", doc = "The maximum reference-genome distance between consecutive heterozygous sites to permit merging phased VCF records into a MNP record", required = false)
    protected int maxGenomicDistanceForMNP = 1;

    @Hidden
    @Argument(fullName = "outputMultipleBaseCountsFile", shortName = "outputMultipleBaseCountsFile", doc = "File to output cases where a single read has multiple bases at the same position [For DEBUGGING purposes only - DO NOT USE!]", required = false)
    protected File outputMultipleBaseCountsFile = null;
    private MultipleBaseCountsWriter outputMultipleBaseCountsWriter = null;

    public void initialize() {
        if (maxPhaseSites <= 2)
            maxPhaseSites = 2; // by definition, must phase a site relative to previous site [thus, 2 in total]

        /*
         Since we cap each base quality (BQ) by its read's mapping quality (MQ) [in Read.updateBaseAndQuality()], then:
         if minBQ > minMQ, then we require that MQ be >= minBQ as well.
         [Otherwise, we end up capping BQ by MQ only AFTER we tried removing bases with BQ < minBQ, which is WRONG!]

         To do this properly, we set: minMQ = max(minMQ, minBQ)
         */
        MIN_MAPPING_QUALITY_SCORE = Math.max(MIN_MAPPING_QUALITY_SCORE, MIN_BASE_QUALITY_SCORE);

        unphasedSiteQueue = new LinkedList<VariantAndReads>();
        partiallyPhasedSites = new CloneableIteratorLinkedList<UnfinishedVariantAndReads>();

        initializeVcfWriter();

        if (variantStatsFilePrefix != null)
            statsWriter = new PhasingQualityStatsWriter(variantStatsFilePrefix);

        if (outputMultipleBaseCountsFile != null)
            outputMultipleBaseCountsWriter = new MultipleBaseCountsWriter(outputMultipleBaseCountsFile);
    }

    private void initializeVcfWriter() {
        // Wrapper VCFWriters will take ownership of inner writers iff: inner writer != origWriter [which wasn't created here]
        VCFWriter origWriter = writer;

        if (enableMergePhasedSegregatingPolymorphismsToMNP)
            writer = new MergeSegregatingAlternateAllelesVCFWriter(writer, getToolkit().getGenomeLocParser(), getToolkit().getArguments().referenceFile, maxGenomicDistanceForMNP, logger, writer != origWriter);

        /* Due to discardIrrelevantPhasedSites(), the startDistance spanned by [partiallyPhasedSites.peek(), unphasedSiteQueue.peek()] is <= cacheWindow
           Due to processQueue(), the startDistance spanned by [unphasedSiteQueue.peek(), mostDownstreamLocusReached] is <= cacheWindow
           Hence, the startDistance between: partiallyPhasedSites.peek() --> mostDownstreamLocusReached is <= 2 * cacheWindow

           Therefore, can write the filtered records located at mostDownstreamLocusReached (if any) to SortingVCFWriter, even though partiallyPhasedSites.peek() has not yet been written.

           But, NOTE that map() is careful to pass out a list of records to be written that FIRST includes any records discarded due to having reached mostDownstreamLocusReached,
           and only THEN records located at mostDownstreamLocusReached.  The opposite order in map() would violate the startDistance limits imposed when contracting SortingVCFWriter with (2 * cacheWindow).
         */
        writer = new SortingVCFWriter(writer, 2 * cacheWindow, writer != origWriter);

        // setup the header fields:
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        // Phasing-specific INFO fields:
        hInfo.add(new VCFFormatHeaderLine(PQ_KEY, 1, VCFHeaderLineType.Float, "Read-backed phasing quality"));
        hInfo.add(new VCFInfoHeaderLine(PHASING_INCONSISTENT_KEY, 0, VCFHeaderLineType.Flag, "Are the reads significantly haplotype-inconsistent?"));

        // todo -- fix samplesToPhase
        String trackName = variantCollection.variants.getName();
        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), Arrays.asList(trackName));
        Set<String> samples = new TreeSet<String>(samplesToPhase == null ? rodNameToHeader.get(trackName).getGenotypeSamples() : samplesToPhase);
        writer.writeHeader(new VCFHeader(hInfo, samples));
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

        mostDownstreamLocusReached = ref.getLocus();
        if (DEBUG) logger.debug("map() at: " + mostDownstreamLocusReached);

        PhasingStats phaseStats = new PhasingStats();
        List<VariantContext> unprocessedList = new LinkedList<VariantContext>();

        for (VariantContext vc : tracker.getValues(variantCollection.variants, context.getLocation())) {
            if (samplesToPhase != null) vc = reduceVCToSamples(vc, samplesToPhase);

            if (ReadBackedPhasingWalker.processVariantInPhasing(vc)) {
                VariantAndReads vr = new VariantAndReads(vc, context);
                unphasedSiteQueue.add(vr);

                if (DEBUG)
                    logger.debug("Added variant to queue = " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant));
            }
            else {
                unprocessedList.add(vc); // Finished with the unprocessed variant, and writer can enforce sorting on-the-fly

                if (DEBUG)
                    logger.debug("Unprocessed variant = " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
            }

            int numReads = 0;
            if (context.hasBasePileup()) {
                numReads = context.getBasePileup().getNumberOfElements();
            }
            else if (context.hasExtendedEventPileup()) {
                numReads = context.getExtendedEventPileup().getNumberOfElements();
            }
            PhasingStats addInPhaseStats = new PhasingStats(numReads, 1);
            phaseStats.addIn(addInPhaseStats);
        }

        List<VariantContext> completedList = processQueue(phaseStats, false);
        completedList.addAll(unprocessedList); // add unprocessedList on to the END of completedList so that the processQueue() results, which are necessarily more upstream, are first!

        return new PhasingStatsAndOutput(phaseStats, completedList);
    }

    private static final Set<String> KEYS_TO_KEEP_IN_REDUCED_VCF = new HashSet<String>(Arrays.asList(PQ_KEY));

    private VariantContext reduceVCToSamples(VariantContext vc, List<String> samplesToPhase) {
//        for ( String sample : samplesToPhase )
//            logger.debug(String.format("  Sample %s has genotype %s, het = %s", sample, vc.getGenotype(sample), vc.getGenotype(sample).isHet() ));
        VariantContext subvc = vc.subContextFromGenotypes(vc.getGenotypes(samplesToPhase).values());
//        logger.debug("original VC = " + vc);
//        logger.debug("sub      VC = " + subvc);
        return VariantContextUtils.pruneVariantContext(subvc, KEYS_TO_KEEP_IN_REDUCED_VCF);
    }

    private List<VariantContext> processQueue(PhasingStats phaseStats, boolean processAll) {
        List<VariantContext> oldPhasedList = new LinkedList<VariantContext>();

        while (!unphasedSiteQueue.isEmpty()) {
            if (!processAll) { // otherwise, phase until the end of unphasedSiteQueue
                VariantContext nextToPhaseVc = unphasedSiteQueue.peek().variant;
                if (startDistancesAreInWindowRange(mostDownstreamLocusReached, VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), nextToPhaseVc))) {
                    /* mostDownstreamLocusReached is still not far enough ahead of nextToPhaseVc to have all phasing information for nextToPhaseVc
                     (note that we ASSUME that the VCF is ordered by <contig,locus>).
                      Note that this will always leave at least one entry (the last one), since mostDownstreamLocusReached is in range of itself.
                    */
                    break;
                }
                // Already saw all variant positions within cacheWindow startDistance ahead of vc (on its contig)
            }
            // Update partiallyPhasedSites before it's used in phaseSite:
            oldPhasedList.addAll(discardIrrelevantPhasedSites());
            if (DEBUG) logger.debug("oldPhasedList(1st) = " + toStringVCL(oldPhasedList));

            VariantAndReads vr = unphasedSiteQueue.remove();
            if (DEBUG)
                logger.debug("Performing phasing for " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant));
            phaseSite(vr, phaseStats);
        }

        // Update partiallyPhasedSites after phaseSite is done:
        oldPhasedList.addAll(discardIrrelevantPhasedSites());
        if (DEBUG) logger.debug("oldPhasedList(2nd) = " + toStringVCL(oldPhasedList));

        if (outputMultipleBaseCountsWriter != null)
            outputMultipleBaseCountsWriter.outputMultipleBaseCounts();

        return oldPhasedList;
    }

    private List<VariantContext> discardIrrelevantPhasedSites() {
        List<VariantContext> vcList = new LinkedList<VariantContext>();

        GenomeLoc nextToPhaseLoc = null;
        if (!unphasedSiteQueue.isEmpty())
            nextToPhaseLoc = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), unphasedSiteQueue.peek().variant);

        while (!partiallyPhasedSites.isEmpty()) {
            if (nextToPhaseLoc != null) { // otherwise, unphasedSiteQueue.isEmpty(), and therefore no need to keep any of the "past"
                UnfinishedVariantAndReads partPhasedVr = partiallyPhasedSites.peek();

                if (startDistancesAreInWindowRange(partPhasedVr.unfinishedVariant.getLocation(), nextToPhaseLoc))
                    // nextToPhaseLoc is still not far enough ahead of partPhasedVr to exclude partPhasedVr from calculations
                    break;
            }
            UnfinishedVariantAndReads uvr = partiallyPhasedSites.remove();
            vcList.add(uvr.unfinishedVariant.toVariantContext());
        }

        return vcList;
    }

    /* Phase vc (removed head of unphasedSiteQueue) using all VariantContext objects in
       partiallyPhasedSites, and all in unphasedSiteQueue that are within cacheWindow startDistance ahead of vc (on its contig).

       ASSUMES: All VariantContexts in unphasedSiteQueue are in positions downstream of vc (head of queue).
     */

    private void phaseSite(VariantAndReads vr, PhasingStats phaseStats) {
        VariantContext vc = vr.variant;
        logger.debug("Will phase vc = " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));

        UnfinishedVariantAndReads uvr = new UnfinishedVariantAndReads(vr);
        UnfinishedVariantContext uvc = uvr.unfinishedVariant;

        // Perform per-sample phasing:
        Map<String, Genotype> sampGenotypes = vc.getGenotypes();
        Map<String, PhaseCounts> samplePhaseStats = new TreeMap<String, PhaseCounts>();
        for (Map.Entry<String, Genotype> sampGtEntry : sampGenotypes.entrySet()) {
            String samp = sampGtEntry.getKey();
            Genotype gt = sampGtEntry.getValue();

            if (DEBUG) logger.debug("sample = " + samp);
            if (isUnfilteredCalledDiploidGenotype(gt)) {
                if (gt.isHom()) { // Note that this Genotype may be replaced later to contain the PQ of a downstream het site that was phased relative to a het site lying upstream of this hom site:
                    // true <-> can trivially phase a hom site relative to ANY previous site:
                    Genotype phasedGt = new Genotype(gt.getSampleName(), gt.getAlleles(), gt.getNegLog10PError(), gt.getFilters(), gt.getAttributes(), true);
                    uvc.setGenotype(samp, phasedGt);
                }
                else if (gt.isHet()) { // Attempt to phase this het genotype relative to the previous het genotype
                    PhasingWindow phaseWindow = new PhasingWindow(vr, samp);
                    if (phaseWindow.hasPreviousHets()) { // Otherwise, nothing to phase this against
                        SNPallelePair allelePair = new SNPallelePair(gt);
                        if (DEBUG) logger.debug("Want to phase TOP vs. BOTTOM for: " + "\n" + allelePair);

                        CloneableIteratorLinkedList.CloneableIterator<UnfinishedVariantAndReads> prevHetAndInteriorIt = phaseWindow.prevHetAndInteriorIt;
                        /* Notes:
                        1. Call to next() advances iterator to next position in partiallyPhasedSites.
                        2. prevHetGenotype != null, since otherwise prevHetAndInteriorIt would not have been chosen to point to its UnfinishedVariantAndReads.
                        */
                        UnfinishedVariantContext prevUvc = prevHetAndInteriorIt.next().unfinishedVariant;
                        Genotype prevHetGenotype = prevUvc.getGenotype(samp);

                        PhaseResult pr = phaseSampleAtSite(phaseWindow);
                        boolean genotypesArePhased = passesPhasingThreshold(pr.phaseQuality);

                        if (pr.phasingContainsInconsistencies) {
                            if (DEBUG)
                                logger.debug("MORE than " + (MAX_FRACTION_OF_INCONSISTENT_READS * 100) + "% of the reads are inconsistent for phasing of " + VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
                            uvc.setPhasingInconsistent();
                        }

                        if (genotypesArePhased) {
                            SNPallelePair prevAllelePair = new SNPallelePair(prevHetGenotype);

                            if (DEBUG)
                                logger.debug("THE PHASE PREVIOUSLY CHOSEN FOR PREVIOUS:\n" + prevAllelePair + "\n");
                            if (DEBUG) logger.debug("THE PHASE CHOSEN HERE:\n" + allelePair + "\n\n");

                            ensurePhasing(allelePair, prevAllelePair, pr.haplotype);
                            Map<String, Object> gtAttribs = new HashMap<String, Object>(gt.getAttributes());
                            gtAttribs.put(PQ_KEY, pr.phaseQuality);
                            Genotype phasedGt = new Genotype(gt.getSampleName(), allelePair.getAllelesAsList(), gt.getNegLog10PError(), gt.getFilters(), gtAttribs, genotypesArePhased);
                            uvc.setGenotype(samp, phasedGt);
                        }

                        // Now, update the 0 or more "interior" hom sites in between the previous het site and this het site:
                        while (prevHetAndInteriorIt.hasNext()) {
                            UnfinishedVariantContext interiorUvc = prevHetAndInteriorIt.next().unfinishedVariant;
                            Genotype handledGt = interiorUvc.getGenotype(samp);
                            if (handledGt == null || !isUnfilteredCalledDiploidGenotype(handledGt))
                                throw new ReviewedStingException("LOGICAL error: should not have breaks WITHIN haplotype");
                            if (!handledGt.isHom())
                                throw new ReviewedStingException("LOGICAL error: should not have anything besides hom sites IN BETWEEN two het sites");

                            // Use the same phasing consistency and PQ for each hom site in the "interior" as for the het-het phase:
                            if (pr.phasingContainsInconsistencies)
                                interiorUvc.setPhasingInconsistent();

                            if (genotypesArePhased) {
                                Map<String, Object> handledGtAttribs = new HashMap<String, Object>(handledGt.getAttributes());
                                handledGtAttribs.put(PQ_KEY, pr.phaseQuality);
                                Genotype phasedHomGt = new Genotype(handledGt.getSampleName(), handledGt.getAlleles(), handledGt.getNegLog10PError(), handledGt.getFilters(), handledGtAttribs, genotypesArePhased);
                                interiorUvc.setGenotype(samp, phasedHomGt);
                            }
                        }

                        if (statsWriter != null)
                            statsWriter.addStat(samp, VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc), startDistance(prevUvc, vc), pr.phaseQuality, phaseWindow.readsAtHetSites.size(), phaseWindow.hetGenotypes.length);

                        PhaseCounts sampPhaseCounts = samplePhaseStats.get(samp);
                        if (sampPhaseCounts == null) {
                            sampPhaseCounts = new PhaseCounts();
                            samplePhaseStats.put(samp, sampPhaseCounts);
                        }
                        sampPhaseCounts.numTestedSites++;

                        if (pr.phasingContainsInconsistencies) {
                            if (genotypesArePhased)
                                sampPhaseCounts.numInconsistentSitesPhased++;
                            else
                                sampPhaseCounts.numInconsistentSitesNotPhased++;
                        }

                        if (genotypesArePhased)
                            sampPhaseCounts.numPhased++;
                    }
                }
            }
        }

        partiallyPhasedSites.add(uvr); // only add it in now, since don't want it to be there during phasing
        phaseStats.addIn(new PhasingStats(samplePhaseStats));
    }

    public boolean passesPhasingThreshold(double PQ) {
        return PQ >= phaseQualityThresh;
    }

    private static class GenotypeAndReadBases {
        public Genotype genotype;
        public ReadBasesAtPosition readBases;
        public GenomeLoc loc;

        public GenotypeAndReadBases(Genotype genotype, ReadBasesAtPosition readBases, GenomeLoc loc) {
            this.genotype = genotype;
            this.readBases = readBases;
            this.loc = loc;
        }
    }

    private class PhasingWindow {
        private Genotype[] hetGenotypes = null;
        private CloneableIteratorLinkedList.CloneableIterator<UnfinishedVariantAndReads> prevHetAndInteriorIt = null;
        private int phasingSiteIndex = -1;
        private Map<String, PhasingRead> readsAtHetSites = null;

        public boolean hasPreviousHets() {
            return phasingSiteIndex > 0;
        }

        // ASSUMES that: isUnfilteredCalledDiploidGenotype(vrGt) && vrGt.isHet() [vrGt = vr.variant.getGenotype(sample)]

        public PhasingWindow(VariantAndReads vr, String sample) {
            List<GenotypeAndReadBases> listHetGenotypes = new LinkedList<GenotypeAndReadBases>();

            // Include previously phased sites in the phasing computation:
            CloneableIteratorLinkedList.CloneableIterator<UnfinishedVariantAndReads> phasedIt = partiallyPhasedSites.iterator();
            while (phasedIt.hasNext()) {
                UnfinishedVariantAndReads phasedVr = phasedIt.next();
                Genotype gt = phasedVr.unfinishedVariant.getGenotype(sample);
                if (gt == null || !isUnfilteredCalledDiploidGenotype(gt)) { // constructed haplotype must start AFTER this "break"
                    listHetGenotypes.clear(); // clear out any history
                }
                else if (gt.isHet()) {
                    GenotypeAndReadBases grb = new GenotypeAndReadBases(gt, phasedVr.sampleReadBases.get(sample), phasedVr.unfinishedVariant.getLocation());
                    listHetGenotypes.add(grb);
                    if (DEBUG) logger.debug("Using UPSTREAM het site = " + grb.loc);
                    prevHetAndInteriorIt = phasedIt.clone();
                }
            }
            phasingSiteIndex = listHetGenotypes.size();
            if (phasingSiteIndex == 0) { // no previous sites against which to phase
                hetGenotypes = null;
                prevHetAndInteriorIt = null;
                return;
            }
            prevHetAndInteriorIt.previous(); // so that it points to the previous het site [and NOT one after it, due to the last call to next()]

            // Add the (het) position to be phased:
            GenomeLoc phaseLocus = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vr.variant);
            GenotypeAndReadBases grbPhase = new GenotypeAndReadBases(vr.variant.getGenotype(sample), vr.sampleReadBases.get(sample), phaseLocus);
            listHetGenotypes.add(grbPhase);
            if (DEBUG)
                logger.debug("PHASING het site = " + grbPhase.loc + " [phasingSiteIndex = " + phasingSiteIndex + "]");

            // Include as-of-yet unphased sites in the phasing computation:
            for (VariantAndReads nextVr : unphasedSiteQueue) {
                if (!startDistancesAreInWindowRange(vr.variant, nextVr.variant)) //nextVr too far ahead of the range used for phasing vc
                    break;
                Genotype gt = nextVr.variant.getGenotype(sample);
                if (gt == null || !isUnfilteredCalledDiploidGenotype(gt)) { // constructed haplotype must end BEFORE this "break"
                    break;
                }
                else if (gt.isHet()) {
                    GenotypeAndReadBases grb = new GenotypeAndReadBases(gt, nextVr.sampleReadBases.get(sample), VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), nextVr.variant));
                    listHetGenotypes.add(grb);
                    if (DEBUG) logger.debug("Using DOWNSTREAM het site = " + grb.loc);
                }
            }

            // First, assemble the "sub-reads" from the COMPLETE WINDOW-BASED SET of heterozygous positions for this sample:
            buildReadsAtHetSites(listHetGenotypes, sample, grbPhase.loc);

            // Remove extraneous reads (those that do not "connect" the two core phasing sites):
            Set<String> onlyKeepReads = removeExtraneousReads(listHetGenotypes.size());

            // Dynamically modify the window to only include sites which have a non-empty set of reads:
            listHetGenotypes = removeExtraneousSites(listHetGenotypes);

            // In any case, must still trim the window size to be "feasible"
            // [**NOTE**: May want to do this to try maximize the preservation of paths from (phasingSiteIndex - 1) to phasingSiteIndex]:
            if (listHetGenotypes.size() > maxPhaseSites) {
                listHetGenotypes = trimWindow(listHetGenotypes, sample, phaseLocus);

                // Can now remove any extra reads (and then sites):
                buildReadsAtHetSites(listHetGenotypes, onlyKeepReads);
                onlyKeepReads = removeExtraneousReads(listHetGenotypes.size());
                listHetGenotypes = removeExtraneousSites(listHetGenotypes);
            }

            // Lastly, assemble the "sub-reads" from the FINAL SET of heterozygous positions for this sample:
            buildReadsAtHetSites(listHetGenotypes, onlyKeepReads);

            // Copy to a fixed-size array:
            if (DEBUG)
                logger.debug("FINAL phasing window of " + listHetGenotypes.size() + " sites:\n" + toStringGRL(listHetGenotypes));
            hetGenotypes = new Genotype[listHetGenotypes.size()];
            int index = 0;
            for (GenotypeAndReadBases copyGrb : listHetGenotypes)
                hetGenotypes[index++] = copyGrb.genotype;
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phasingLoc) {
            buildReadsAtHetSites(listHetGenotypes, sample, phasingLoc, null);
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, Set<String> onlyKeepReads) {
            buildReadsAtHetSites(listHetGenotypes, null, null, onlyKeepReads);
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phasingLoc, Set<String> onlyKeepReads) {
            readsAtHetSites = new HashMap<String, PhasingRead>();

            int index = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                ReadBasesAtPosition readBases = grb.readBases;
                if (readBases != null) {
                    for (ReadBase rb : readBases) {
                        String readName = rb.readName;
                        if (onlyKeepReads != null && !onlyKeepReads.contains(readName)) // if onlyKeepReads exists, ignore reads not in onlyKeepReads
                            continue;

                        PhasingRead rd = readsAtHetSites.get(readName);
                        if (rd == null) {
                            rd = new PhasingRead(listHetGenotypes.size(), rb.mappingQual);
                            readsAtHetSites.put(readName, rd);
                        }
                        else if (outputMultipleBaseCountsWriter != null && rd.getBase(index) != null // rd already has a base at index
                                && sample != null && phasingLoc != null) {
                            outputMultipleBaseCountsWriter.setMultipleBases(new SampleReadLocus(sample, readName, grb.loc), phasingLoc, rd.getBase(index), rb.base);
                        }

                        // Arbitrarily updates to the last base observed for this sample and read (rb.base):
                        rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
                    }
                }
                index++;
            }
            if (DEBUG) logger.debug("Number of sites in window = " + index);

            if (DEBUG && logger.isDebugEnabled()) {
                logger.debug("ALL READS [phasingSiteIndex = " + phasingSiteIndex + "]:");
                for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                    String rdName = nameToReads.getKey();
                    PhasingRead rd = nameToReads.getValue();
                    logger.debug(rd + "\t" + rdName);
                }
            }
        }

        private class EdgeToReads {
            private TreeMap<PhasingGraphEdge, List<String>> edgeReads;

            public EdgeToReads() {
                this.edgeReads = new TreeMap<PhasingGraphEdge, List<String>>(); // implemented GraphEdge.compareTo()
            }

            public void addRead(PhasingGraphEdge e, String readName) {
                List<String> reads = edgeReads.get(e);
                if (reads == null) {
                    reads = new LinkedList<String>();
                    edgeReads.put(e, reads);
                }
                reads.add(readName);
            }

            public List<String> getReads(PhasingGraphEdge e) {
                return edgeReads.get(e);
            }
        }

        private class IntegerSet implements Iterable<Integer> {
            private Set<Integer> list;

            public IntegerSet(Set<Integer> list) {
                this.list = list;
            }

            public boolean contains(int i) {
                return list.contains(i);
            }

            public Iterator<Integer> iterator() {
                return list.iterator();
            }

            public String toString() {
                StringBuilder sb = new StringBuilder();
                for (int i : this) {
                    sb.append(i + ", ");
                }
                return sb.toString();
            }
        }

        public Set<String> removeExtraneousReads(int numHetSites) {
            PhasingGraph readGraph = new PhasingGraph(numHetSites);
            EdgeToReads edgeToReads = new EdgeToReads();
            Set<Integer> sitesWithEdges = new TreeSet<Integer>();

            for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                PhasingRead rd = nameToReads.getValue();

                int[] siteInds = rd.getNonNullIndices();
                // Connect each pair of non-null sites in rd:
                for (int i = 0; i < siteInds.length; i++) {
                    for (int j = i + 1; j < siteInds.length; j++) {
                        PhasingGraphEdge e = new PhasingGraphEdge(siteInds[i], siteInds[j]);
                        if (DEBUG) logger.debug("Read = " + rdName + " is adding edge: " + e);
                        readGraph.addEdge(e);

                        edgeToReads.addRead(e, rdName);

                        sitesWithEdges.add(e.getV1());
                        sitesWithEdges.add(e.getV2());
                    }
                }
            }
            if (DEBUG) logger.debug("Read graph:\n" + readGraph);
            Set<String> keepReads = new HashSet<String>();

            /* Check which Reads are involved in acyclic paths from (phasingSiteIndex - 1) to (phasingSiteIndex):

               In detail:
               Every Read links EACH pair of sites for which it contains bases.  Then, each such edge is added to a "site connectivity graph".
               A read provides non-trivial bias toward the final haplotype decision if it participates in a path from prev ---> cur.  This is tested by
               considering each edge that the read contributes.  For edge e=(v1,v2), if there exists a path from prev ---> v1 [that doesn't include v2] and
               cur ---> v2 [that doesn't include v1], then there is a path from prev ---> cur that uses e, hence making the read significant.
               By excluding each vertex's edges and then calculating connected components, we are able to make the determination, for example,
               if a path exists from prev ---> v1 that excludes v2.

               Furthermore, if the path DOES use other edges that exist solely due to the read, then that's fine, since adding in the read will give those edges as well.
               And, if the path uses edges from other reads, then keeping all other reads that contribute those edges
               [which will happen since those edges are also in paths from prev ---> cur] is sufficient for this path to exist.

               NOTE:
               If we would use NON-UNIFORM priors for the various haplotypes consistent with a margnialized haplotype, then this calculation would not be correct, since the equivalence of:
               1. The read affects the final marginal haplotype posterior probability (for general mapping and base quality values).
               2. The read has edges involved in a path from prev ---> cur.
               DEPENDS STRONGLY on the fact that all haplotypes have the same EXACT prior.

               This is due to the following:
               [We denote:
               R = set of all reads
               r = a single read
               "AA + CC" = AA on top chromosome, CC on bottom chromosome]

               Note that since there are only two haplotype possibilities:
               P(AA + CC | R) + P(AC + CA | R) = 1

               Now, if we assume that all haplotypes consistent with AA + CC have the same prior probability [P(AA + CC | R)], then:
               P(AA + CC | R)
               = P(AAAA + CCCC | R) + ... + P(AACC + CCAA | R)
               = [P(AAAA + CCCC , R) + ... + P(AACC + CCAA , R)] / P(R)
               \propto P(AAAA + CCCC , R) + ... + P(AACC + CCAA , R)
               = P(R | AAAA + CCCC)*P(AAAA + CCCC) + ... + P(R | AACC + CCAA)*P(AACC + CCAA)
               = P(AA + CC | R) * [P(R | AAAA + CCCC) + ... + P(R | AACC + CCAA)]
               
               Since we assume independence between reads given a particular haplotype [P(R | AAAA + CCCC) = \prod_r P(r | AAAA + CCCC)],
               a new read r affects P(AA + CC | R) by multiplying each of the terms in the sum by, e.g., P(r | AAAA + CCCC).
               Therefore, if these values do not affect the ratio of:
               (I) [P(R | AAAA + CCCC) + ... + P(R | AACC + CCAA)] / [P(R | ACAA + CACC) + ... + P(R | ACCC + CAAA)]
               then they do not affect the value of:
               (II) P(AA + CC | R) / P(AC + CA | R)   [which uniquely defines their values, since they sum to 1]

               And, the P(r | AAAA + CCCC), ..., P(r | ACCC + CAAA) do not affect ratio (I) iff r's edges do not take part in a path from prev to cur in combination with the other reads in R.
             */
            int prev = phasingSiteIndex - 1;
            int cur = phasingSiteIndex;

            if (!readGraph.getConnectedComponents().inSameSet(prev, cur)) { // There is NO path between cur and prev
                if (DEBUG)
                    logger.debug("NO READ PATH between PHASE site [" + cur + "] and UPSTREAM site [" + prev + "]");
                readsAtHetSites.clear();
                return keepReads;
            }

            /* Check the connected components of prev and cur when removing each individual vertex's edges:
               [Total run-time: for each vertex, calculate connected components after removing it's edges: O(V * E)]
             */
            IntegerSet[] removedSiteSameCCAsPrev = new IntegerSet[numHetSites];
            IntegerSet[] removedSiteSameCCAsCur = new IntegerSet[numHetSites];
            for (int i : sitesWithEdges) {
                if (DEBUG) logger.debug("Calculating CC after removing edges of site: " + i);

                // Remove all edges incident to i and see which positions have paths to prev and cur:
                Collection<PhasingGraphEdge> removedEdges = readGraph.removeAllIncidentEdges(i);

                // Run-time for efficiently calculating connected components using DisjointSet: O(E)
                DisjointSet ccAfterRemove = readGraph.getConnectedComponents();
                removedSiteSameCCAsPrev[i] = new IntegerSet(ccAfterRemove.inSameSetAs(prev, sitesWithEdges));
                removedSiteSameCCAsCur[i] = new IntegerSet(ccAfterRemove.inSameSetAs(cur, sitesWithEdges));

                if (DEBUG) logger.debug("Same CC as previous [" + prev + "]: " + removedSiteSameCCAsPrev[i]);
                if (DEBUG) logger.debug("Same CC as current  [" + cur + "]: " + removedSiteSameCCAsCur[i]);

                // Add the removed edges back in:
                readGraph.addEdges(removedEdges);
            }

            for (PhasingGraphEdge e : readGraph) {
                if (DEBUG) logger.debug("Testing the path-connectivity of Edge: " + e);

                /* Edge e={v1,v2} contributes a path between prev and cur for testRead iff:
                   testRead[v1] != null, testRead[v2] != null, and there is a path from prev ---> v1 -> v2 ---> cur  [or vice versa].
                   Note that the path from prev ---> v1 will NOT contain v2, since we removed all of v2's edges,
                   and the path from v2 ---> cur will NOT contain v1.
                 */
                boolean prevTo2and1ToCur = removedSiteSameCCAsPrev[e.getV1()].contains(e.getV2()) && removedSiteSameCCAsCur[e.getV2()].contains(e.getV1());
                boolean prevTo1and2ToCur = removedSiteSameCCAsPrev[e.getV2()].contains(e.getV1()) && removedSiteSameCCAsCur[e.getV1()].contains(e.getV2());

                if (prevTo2and1ToCur || prevTo1and2ToCur) {
                    for (String readName : edgeToReads.getReads(e)) {
                        keepReads.add(readName);

                        if (DEBUG && logger.isDebugEnabled()) {
                            if (prevTo2and1ToCur)
                                logger.debug("Keep read " + readName + " due to path: " + prev + " ---> " + e.getV2() + " -> " + e.getV1() + " ---> " + cur);
                            else
                                logger.debug("Keep read " + readName + " due to path: " + prev + " ---> " + e.getV1() + " -> " + e.getV2() + " ---> " + cur);
                        }
                    }
                }
            }

            // Retain only the reads that contain an edge in a path connecting prev and cur:
            Iterator<Map.Entry<String, PhasingRead>> readIt = readsAtHetSites.entrySet().iterator();
            while (readIt.hasNext()) {
                Map.Entry<String, PhasingRead> nameToReads = readIt.next();
                String rdName = nameToReads.getKey();
                if (!keepReads.contains(rdName)) {
                    readIt.remove();
                    if (DEBUG) logger.debug("Removing extraneous read: " + rdName);
                }
            }

            return keepReads;
        }

        private List<GenotypeAndReadBases> removeExtraneousSites(List<GenotypeAndReadBases> listHetGenotypes) {
            Set<Integer> sitesWithReads = new HashSet<Integer>();
            for (Map.Entry<String, PhasingRead> nameToReads : readsAtHetSites.entrySet()) {
                PhasingRead rd = nameToReads.getValue();
                for (int i : rd.getNonNullIndices())
                    sitesWithReads.add(i);
            }

            // Remove all sites that have no read bases:
            List<GenotypeAndReadBases> keepHetSites = new LinkedList<GenotypeAndReadBases>();
            int index = 0;
            int numPrecedingRemoved = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                boolean keepSite = sitesWithReads.contains(index);
                if (DEBUG && logger.isDebugEnabled() && !keepSite)
                    logger.debug("Removing read-less site " + grb.loc);

                if (keepSite || index == phasingSiteIndex || index == phasingSiteIndex - 1) {
                    keepHetSites.add(grb);
                    if (!keepSite)
                        if (DEBUG)
                            logger.debug("Although current or previous sites have no relevant reads, continuing empty attempt to phase them [for sake of program flow]...");
                }
                else if (index <= phasingSiteIndex)
                    numPrecedingRemoved++;

                index++;
            }

            phasingSiteIndex -= numPrecedingRemoved;
            return keepHetSites;
        }

        private List<GenotypeAndReadBases> trimWindow(List<GenotypeAndReadBases> listHetGenotypes, String sample, GenomeLoc phaseLocus) {
            if (DEBUG)
                logger.warn("Trying to phase sample " + sample + " at locus " + phaseLocus + " within a window of " + cacheWindow + " bases yields " + listHetGenotypes.size() + " heterozygous sites to phase:\n" + toStringGRL(listHetGenotypes));

            int prevSiteIndex = phasingSiteIndex - 1; // index of previous in listHetGenotypes
            int numToUse = maxPhaseSites - 2; // since always keep previous and current het sites!

            int numOnLeft = prevSiteIndex;
            int numOnRight = listHetGenotypes.size() - (phasingSiteIndex + 1);

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
            listHetGenotypes = listHetGenotypes.subList(startIndex, stopIndex);
            if (DEBUG)
                logger.warn("NAIVELY REDUCED to " + listHetGenotypes.size() + " sites:\n" + toStringGRL(listHetGenotypes));

            return listHetGenotypes;
        }
    }

    private PhaseResult phaseSampleAtSite(PhasingWindow phaseWindow) {
        /* Will map a phase and its "complement" to a single representative phase,
          and marginalizeAsNewTable() marginalizes to 2 positions [starting at the previous position, and then the current position]:
        */
        HaplotypeTableCreator tabCreator = new TableCreatorOfHaplotypeAndComplementForDiploidAlleles(phaseWindow.hetGenotypes, phaseWindow.phasingSiteIndex - 1, 2);
        PhasingTable sampleHaps = tabCreator.getNewTable();

        if (DEBUG && logger.isDebugEnabled()) {
            logger.debug("Number of USED reads [connecting the two positions to be phased] at sites: " + phaseWindow.readsAtHetSites.size());
            logger.debug("USED READS:");
            for (Map.Entry<String, PhasingRead> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                PhasingRead rd = nameToReads.getValue();
                logger.debug(rd + "\t" + rdName);
            }
        }

        // Update the phasing table based on each of the sub-reads for this sample:
        MaxHaplotypeAndQuality prevMaxHapAndQual = null;

        int numHighQualityIterations = 0;
        int numInconsistentIterations = 0;

        double totalAbsPQchange = 0;
        int numPQchangesObserved = 0;

        for (Map.Entry<String, PhasingRead> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
            PhasingRead rd = nameToReads.getValue();
            if (DEBUG) logger.debug("\nrd = " + rd + "\tname = " + nameToReads.getKey());

            for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
                PhasingScore score = rd.matchHaplotypeClassScore(pte.getHaplotypeClass());
                pte.getScore().integrateReadScore(score);
                if (DEBUG) logger.debug("score(" + rd + ", " + pte.getHaplotypeClass() + ") = " + score);
            }

            // Check the current best haplotype assignment and compare it to the previous one:
            MaxHaplotypeAndQuality curMaxHapAndQual = new MaxHaplotypeAndQuality(sampleHaps, false);
            if (DEBUG)
                logger.debug("CUR MAX hap:\t" + curMaxHapAndQual.maxEntry.getHaplotypeClass() + "\tcurPhaseQuality:\t" + curMaxHapAndQual.phaseQuality);
            if (prevMaxHapAndQual != null) {
                double changeInPQ = prevMaxHapAndQual.phaseQuality - curMaxHapAndQual.phaseQuality;

                if (passesPhasingThreshold(prevMaxHapAndQual.phaseQuality)) {
                    numHighQualityIterations++;
                    if (!curMaxHapAndQual.hasSameRepresentativeHaplotype(prevMaxHapAndQual) || // switched phase
                            (numPQchangesObserved > 0 && changeInPQ > FRACTION_OF_MEAN_PQ_CHANGES * (totalAbsPQchange / numPQchangesObserved))) { // a "significant" decrease in PQ
                        if (DEBUG) logger.debug("Inconsistent read found!");
                        numInconsistentIterations++;
                    }
                }

                totalAbsPQchange += Math.abs(changeInPQ);
                numPQchangesObserved++;
            }
            prevMaxHapAndQual = curMaxHapAndQual;
        }

        if (DEBUG) logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");
        MaxHaplotypeAndQuality maxHapQual = new MaxHaplotypeAndQuality(sampleHaps, true);
        double posteriorProb = maxHapQual.maxEntry.getScore().getValue();

        if (DEBUG)
            logger.debug("MAX hap:\t" + maxHapQual.maxEntry.getHaplotypeClass() + "\tposteriorProb:\t" + posteriorProb + "\tphaseQuality:\t" + maxHapQual.phaseQuality);
        if (DEBUG)
            logger.debug("Number of used reads " + phaseWindow.readsAtHetSites.size() + "; number of high PQ iterations " + numHighQualityIterations + "; number of inconsistencies " + numInconsistentIterations);

        boolean phasingContainsInconsistencies = false;
        if (numInconsistentIterations / (double) numHighQualityIterations > MAX_FRACTION_OF_INCONSISTENT_READS)
            phasingContainsInconsistencies = true;

        return new PhaseResult(maxHapQual.getRepresentative(), maxHapQual.phaseQuality, phasingContainsInconsistencies);
    }

    private static class MaxHaplotypeAndQuality {
        public PhasingTable.PhasingTableEntry maxEntry;
        public double phaseQuality;

        public MaxHaplotypeAndQuality(PhasingTable hapTable, boolean printDebug) {
            // Marginalize each haplotype to its first 2 positions:
            hapTable = HaplotypeTableCreator.marginalizeAsNewTable(hapTable);
            if (DEBUG && printDebug)
                logger.debug("\nPhasing table [AFTER MAPPING]:\n" + hapTable + "\n");

            calculateMaxHapAndPhasingQuality(hapTable, printDebug);
        }

        // Calculates maxEntry and its PQ (within table hapTable):

        private void calculateMaxHapAndPhasingQuality(PhasingTable hapTable, boolean printDebug) {
            hapTable.normalizeScores();
            if (DEBUG && printDebug)
                logger.debug("\nPhasing table [AFTER NORMALIZATION]:\n" + hapTable + "\n");

            // Determine the phase at this position:
            this.maxEntry = hapTable.maxEntry();

            // convert posteriorProb to PHRED scale, but do NOT cap the quality as in QualityUtils.probToQual(posteriorProb):
            PreciseNonNegativeDouble sumErrorProbs = new PreciseNonNegativeDouble(ZERO);
            for (PhasingTable.PhasingTableEntry pte : hapTable) {
                if (pte != maxEntry)
                    sumErrorProbs.plusEqual(pte.getScore());
            }
            this.phaseQuality = -10.0 * (sumErrorProbs.getLog10Value());
        }

        public boolean hasSameRepresentativeHaplotype(MaxHaplotypeAndQuality that) {
            return this.getRepresentative().equals(that.getRepresentative());
        }

        private Haplotype getRepresentative() {
            return maxEntry.getHaplotypeClass().getRepresentative();
        }
    }

    /*
        Ensure that curAllelePair is phased relative to prevAllelePair as specified by hap.
     */

    public static void ensurePhasing(SNPallelePair curAllelePair, SNPallelePair prevAllelePair, Haplotype hap) {
        if (hap.size() < 2)
            throw new ReviewedStingException("LOGICAL ERROR: Only considering haplotypes of length > 2!");

        byte prevBase = hap.getBase(0); // The 1st base in the haplotype
        byte curBase = hap.getBase(1);  // The 2nd base in the haplotype

        boolean chosePrevTopChrom = prevAllelePair.matchesTopBase(prevBase);
        boolean choseCurTopChrom = curAllelePair.matchesTopBase(curBase);
        if (chosePrevTopChrom != choseCurTopChrom)
            curAllelePair.swapAlleles();
    }

    private boolean startDistancesAreInWindowRange(VariantContext vc1, VariantContext vc2) {
        return startDistancesAreInWindowRange(VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc1), VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc2));
    }

    private boolean startDistancesAreInWindowRange(GenomeLoc loc1, GenomeLoc loc2) {
        return loc1.distance(loc2) <= cacheWindow; // distance() checks: loc1.onSameContig(loc2)
    }

    private int startDistance(UnfinishedVariantContext uvc1, VariantContext vc2) {
        return uvc1.getLocation().distance(VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc2));
    }

    public PhasingStats reduce(PhasingStatsAndOutput statsAndList, PhasingStats stats) {
        if (statsAndList != null) {
            writeVcList(statsAndList.output);
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
        List<VariantContext> finalList = processQueue(result, true); // process all remaining data
        writeVcList(finalList);
        writer.close();

        if (statsWriter != null)
            statsWriter.close();

        if (outputMultipleBaseCountsWriter != null)
            outputMultipleBaseCountsWriter.close();

        System.out.println("Coverage over ALL samples:");
        System.out.println("Number of reads observed: " + result.getNumReads());
        System.out.println("Number of variant sites observed: " + result.getNumVarSites());
        System.out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));

        System.out.println("\n--- Phasing summary [minimal haplotype quality (PQ): " + phaseQualityThresh + ", maxPhaseSites: " + maxPhaseSites + ", cacheWindow: " + cacheWindow + "] ---");
        for (Map.Entry<String, PhaseCounts> sampPhaseCountEntry : result.getPhaseCounts()) {
            PhaseCounts pc = sampPhaseCountEntry.getValue();
            System.out.print("Sample: " + sampPhaseCountEntry.getKey() + "\tSites tested: " + pc.numTestedSites + "\tSites phased: " + pc.numPhased);
            System.out.println("\tPhase-inconsistent sites: " + (pc.numInconsistentSitesPhased + pc.numInconsistentSitesNotPhased) + " [phased: " + pc.numInconsistentSitesPhased + ", unphased:" + pc.numInconsistentSitesNotPhased + "]");
        }
        System.out.println("");
    }

    private void writeVcList(List<VariantContext> varContList) {
        for (VariantContext vc : varContList)
            writeVCF(vc);
    }

    private void writeVCF(VariantContext vc) {
        if (samplesToPhase == null || vc.isNotFiltered())
            //if ( samplesToPhase == null || (vc.isVariant() && vc.isNotFiltered())) // if we are only operating on specific samples, don't write out all sites, just those where the VC is variant
            WriteVCF.writeVCF(vc, writer, logger);
    }

    public static boolean processVariantInPhasing(VariantContext vc) {
        return vc.isNotFiltered() && ((vc.isSNP() && vc.isBiallelic()) || !vc.isVariant()); // we can handle the non-variant case as well
        //return isUnfilteredBiallelicSNP(vc);
    }


    /*
      Inner classes:
    */

    private class VariantAndReads {
        public VariantContext variant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;

        public VariantAndReads(VariantContext variant, HashMap<String, ReadBasesAtPosition> sampleReadBases) {
            this.variant = variant;
            this.sampleReadBases = sampleReadBases;
        }

        public VariantAndReads(VariantContext variant, AlignmentContext alignment) {
            this.variant = variant;
            this.sampleReadBases = new HashMap<String, ReadBasesAtPosition>();

            if (alignment != null) {
                ReadBackedPileup pileup = null;
                if (alignment.hasBasePileup()) {
                    pileup = alignment.getBasePileup();
                }
                else if (alignment.hasExtendedEventPileup()) {
                    pileup = alignment.getExtendedEventPileup();
                }
                if (pileup != null) {
                    // filter the read-base pileup based on min base and mapping qualities:
                    pileup = pileup.getBaseAndMappingFilteredPileup(MIN_BASE_QUALITY_SCORE, MIN_MAPPING_QUALITY_SCORE);
                    if (pileup != null) {
                        for (final String sample : pileup.getSamples()) {
                            ReadBackedPileup samplePileup = pileup.getPileupForSample(sample);
                            ReadBasesAtPosition readBases = new ReadBasesAtPosition();
                            for (PileupElement p : samplePileup) {
                                if (!p.isDeletion()) // IGNORE deletions for now
                                    readBases.putReadBase(p);
                            }
                            sampleReadBases.put(sample, readBases);
                        }
                    }
                }
            }
        }
    }

    private class UnfinishedVariantAndReads {
        public UnfinishedVariantContext unfinishedVariant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;

        public UnfinishedVariantAndReads(VariantAndReads vr) {
            this.unfinishedVariant = new UnfinishedVariantContext(vr.variant);
            this.sampleReadBases = vr.sampleReadBases;
        }
    }

    // COULD replace with MutableVariantContext if it worked [didn't throw exceptions when trying to call its set() methods]...

    private class UnfinishedVariantContext implements HasGenomeLocation {
        private String name;
        private String contig;
        private int start;
        private int stop;
        private Collection<Allele> alleles;
        private Map<String, Genotype> genotypes;
        private double negLog10PError;
        private Set<String> filters;
        private Map<String, Object> attributes;

        public UnfinishedVariantContext(VariantContext vc) {
            this.name = vc.getSource();
            this.contig = vc.getChr();
            this.start = vc.getStart();
            this.stop = vc.getEnd();
            this.alleles = vc.getAlleles();
            this.genotypes = new HashMap<String, Genotype>(vc.getGenotypes()); // since vc.getGenotypes() is unmodifiable
            this.negLog10PError = vc.getNegLog10PError();
            this.filters = vc.filtersWereApplied() ? vc.getFilters() : null;
            this.attributes = new HashMap<String, Object>(vc.getAttributes());
        }

        public VariantContext toVariantContext() {
            return new VariantContext(name, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes);
        }

        public GenomeLoc getLocation() {
            return getToolkit().getGenomeLocParser().createGenomeLoc(contig, start, stop);
        }

        public Genotype getGenotype(String sample) {
            return genotypes.get(sample);
        }

        public void setGenotype(String sample, Genotype newGt) {
            genotypes.put(sample, newGt);
        }

        public void setPhasingInconsistent() {
            attributes.put(PHASING_INCONSISTENT_KEY, true);
        }
    }

    private static String toStringGRL(List<GenotypeAndReadBases> grbList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (GenotypeAndReadBases grb : grbList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(grb.loc);
        }
        return sb.toString();
    }

    private String toStringVCL(List<VariantContext> vcList) {
        boolean first = true;
        StringBuilder sb = new StringBuilder();
        for (VariantContext vc : vcList) {
            if (first)
                first = false;
            else
                sb.append(" -- ");

            sb.append(VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), vc));
        }
        return sb.toString();
    }

//
// THIS IMPLEMENTATION WILL FAIL WHEN NOT DEALING WITH SNP Alleles (e.g., MNP or INDEL), SINCE THEN THE Allele.getBases()
// FUNCTION WILL RETURN VARIABLE-LENGTH Byte ARRAYS.  IN THAT CASE, BaseArray/Haplotype/Read WILL NEED TO BE REPLACED WITH
// AN ArrayList OF Allele [OR SIMILAR OBJECT], and WON'T USE: getSingleBase(alleleI)
//

    private static abstract class HaplotypeTableCreator {
        protected Genotype[] genotypes;

        public HaplotypeTableCreator(Genotype[] hetGenotypes) {
            this.genotypes = hetGenotypes;
        }

        abstract public PhasingTable getNewTable();

        protected List<Haplotype> getAllHaplotypes() {
            int numSites = genotypes.length;
            int[] genotypeCards = new int[numSites];
            for (int i = 0; i < numSites; i++)
                genotypeCards[i] = genotypes[i].getPloidy();

            LinkedList<Haplotype> allHaps = new LinkedList<Haplotype>();
            CardinalityCounter alleleCounter = new CardinalityCounter(genotypeCards);
            for (int[] alleleInds : alleleCounter) {
                byte[] hapBases = new byte[numSites];
                for (int i = 0; i < numSites; i++) {
                    Allele alleleI = genotypes[i].getAllele(alleleInds[i]);
                    hapBases[i] = SNPallelePair.getSingleBase(alleleI);
                }
                allHaps.add(new Haplotype(hapBases));
            }
            return allHaps;
        }

        public static PhasingTable marginalizeAsNewTable(PhasingTable table) {
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

    private static class TableCreatorOfHaplotypeAndComplementForDiploidAlleles extends HaplotypeTableCreator {
        private SNPallelePair[] SNPallelePairs;
        private int startIndex;
        private int marginalizeLength;

        public TableCreatorOfHaplotypeAndComplementForDiploidAlleles(Genotype[] hetGenotypes, int startIndex, int marginalizeLength) {
            super(hetGenotypes);

            this.SNPallelePairs = new SNPallelePair[genotypes.length];
            for (int i = 0; i < genotypes.length; i++)
                SNPallelePairs[i] = new SNPallelePair(genotypes[i]);

            this.startIndex = startIndex;
            this.marginalizeLength = marginalizeLength;
        }

        public PhasingTable getNewTable() {
            PhasingTable table = new PhasingTable();
            for (Haplotype hap : getAllHaplotypes()) {
                if (SNPallelePairs[startIndex].matchesTopBase(hap.getBase(startIndex))) {
                    /* hap is the "representative" haplotype [DEFINED here to be
                      the one with the top base at the startIndex position.
                      NOTE that it is CRITICAL that this definition be consistent with the representative sub-haplotypes defined below!]
                    */
                    ArrayList<Haplotype> hapList = new ArrayList<Haplotype>();
                    hapList.add(hap);
                    hapList.add(complement(hap));

                    // want marginalizeLength positions starting at startIndex:
                    Haplotype rep = hap.subHaplotype(startIndex, startIndex + marginalizeLength);
                    double hapClassPrior = getHaplotypeRepresentativePrior(rep); // Note that prior is ONLY a function of the representative haplotype

                    HaplotypeClass hapClass = new HaplotypeClass(hapList, rep);
                    table.addEntry(hapClass, hapClassPrior);
                }
            }
            return table;
        }

        // Can change later to weight the representative Haplotypes differently:

        private double getHaplotypeRepresentativePrior(Haplotype rep) {
            return 1.0;
        }

        private Haplotype complement(Haplotype hap) {
            int numSites = SNPallelePairs.length;
            if (hap.size() != numSites)
                throw new ReviewedStingException("INTERNAL ERROR: hap.size() != numSites");

            // Take the other base at EACH position of the Haplotype:
            byte[] complementBases = new byte[numSites];
            for (int i = 0; i < numSites; i++)
                complementBases[i] = SNPallelePairs[i].getOtherBase(hap.getBase(i));

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
        public boolean phasingContainsInconsistencies;

        public PhaseResult(Haplotype haplotype, double phaseQuality, boolean phasingContainsInconsistencies) {
            this.haplotype = haplotype;
            this.phaseQuality = phaseQuality;
            this.phasingContainsInconsistencies = phasingContainsInconsistencies;
        }
    }

    public static boolean isUnfilteredBiallelicSNP(VariantContext vc) {
        return (vc.isNotFiltered() && vc.isSNP() && vc.isBiallelic());
    }

    public static boolean isUnfilteredCalledDiploidGenotype(Genotype gt) {
        return (gt.isNotFiltered() && gt.isCalled() && gt.getPloidy() == 2);
    }

    private class MultipleBaseCountsWriter {
        private BufferedWriter writer = null;
        private TreeMap<SampleReadLocus, MultipleBaseCounts> multipleBaseCounts = null;

        public MultipleBaseCountsWriter(File outputMultipleBaseCountsFile) {
            FileOutputStream output;
            try {
                output = new FileOutputStream(outputMultipleBaseCountsFile);
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Unable to create multiple base count file at location: " + outputMultipleBaseCountsFile);
            }
            this.writer = new BufferedWriter(new OutputStreamWriter(output));

            this.multipleBaseCounts = new TreeMap<SampleReadLocus, MultipleBaseCounts>(); // implemented SampleReadLocus.compareTo()
        }

        public void setMultipleBases(SampleReadLocus srl, GenomeLoc phasingLoc, byte prevBase, byte newBase) {
            MultipleBaseCounts mbc = multipleBaseCounts.get(srl);
            if (mbc == null) {
                mbc = new MultipleBaseCounts(phasingLoc);
                mbc.incrementBaseCount(prevBase); // only now, do we know to note this
                multipleBaseCounts.put(srl, mbc);
            }
            if (mbc.samePhasingLocAs(phasingLoc)) // otherwise, don't want to count these multiple base counts again
                mbc.incrementBaseCount(newBase);

        }

        public void outputMultipleBaseCounts() {
            GenomeLoc nextToPhaseLoc = null;
            if (!unphasedSiteQueue.isEmpty())
                nextToPhaseLoc = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(), unphasedSiteQueue.peek().variant);

            outputMultipleBaseCounts(nextToPhaseLoc);
        }

        private void outputMultipleBaseCounts(GenomeLoc nextToPhaseLoc) {
            try {
                Iterator<Map.Entry<SampleReadLocus, MultipleBaseCounts>> multBaseCountIt = multipleBaseCounts.entrySet().iterator();
                while (multBaseCountIt.hasNext()) {
                    Map.Entry<SampleReadLocus, MultipleBaseCounts> sampleReadLocBaseCountsEntry = multBaseCountIt.next();
                    SampleReadLocus srl = sampleReadLocBaseCountsEntry.getKey();
                    if (nextToPhaseLoc == null || !startDistancesAreInWindowRange(srl.getLocus(), nextToPhaseLoc)) {
                        // Done with entry, so print it and remove it from map:
                        writer.write(srl + "\t" + sampleReadLocBaseCountsEntry.getValue() + "\n");
                        multBaseCountIt.remove();
                    }
                }
                writer.flush();
            } catch (IOException e) {
                throw new RuntimeException("Unable to write to outputMultipleBaseCountsFile", e);
            }
        }

        public void close() {
            outputMultipleBaseCounts(null);

            try {
                writer.flush();
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Unable to close outputMultipleBaseCountsFile");
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
                sb.append(" + ");

            sb.append(h);
        }
        sb.append(" [").append(rep).append("]");
        return sb.toString();
    }
}

class PhasingStats {
    private int numReads;
    private int numVarSites;

    // Map of: sample -> PhaseCounts:
    private Map<String, PhaseCounts> samplePhaseStats;

    public PhasingStats() {
        this(new TreeMap<String, PhaseCounts>());
    }

    public PhasingStats(int numReads, int numVarSites) {
        this.numReads = numReads;
        this.numVarSites = numVarSites;
        this.samplePhaseStats = new TreeMap<String, PhaseCounts>();
    }

    public PhasingStats(Map<String, PhaseCounts> samplePhaseStats) {
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
    public int numInconsistentSitesPhased;
    public int numInconsistentSitesNotPhased;
    public int numPhased;

    public PhaseCounts() {
        this.numTestedSites = 0;
        this.numInconsistentSitesPhased = 0;
        this.numInconsistentSitesNotPhased = 0;
        this.numPhased = 0;
    }

    public void addIn(PhaseCounts other) {
        this.numTestedSites += other.numTestedSites;
        this.numInconsistentSitesPhased += other.numInconsistentSitesPhased;
        this.numInconsistentSitesNotPhased += other.numInconsistentSitesNotPhased;
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

    public void addStat(String sample, GenomeLoc locus, int startDistanceFromPrevious, double phasingQuality, int numReads, int windowSize) {
        BufferedWriter sampWriter = sampleToStatsWriter.get(sample);
        if (sampWriter == null) {
            String fileName = variantStatsFilePrefix + "." + sample + ".locus_distance_PQ_numReads_windowSize.txt";

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
            sampWriter.write(locus + "\t" + startDistanceFromPrevious + "\t" + phasingQuality + "\t" + numReads + "\t" + windowSize + "\n");
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

class SampleReadLocus implements Comparable<SampleReadLocus> {
    private String sample;
    private String read;
    private GenomeLoc locus;

    public SampleReadLocus(String sample, String read, GenomeLoc locus) {
        this.sample = sample;
        this.read = read;
        this.locus = locus;
    }

    public GenomeLoc getLocus() {
        return locus;
    }

    public int compareTo(SampleReadLocus that) {
        int comp = this.sample.compareTo(that.sample);
        if (comp != 0)
            return comp;

        comp = this.read.compareTo(that.read);
        if (comp != 0)
            return comp;

        return this.locus.compareTo(that.locus);
    }

    public String toString() {
        return "Sample " + sample + ", read " + read + ", locus " + locus;
    }
}

class MultipleBaseCounts {
    private Map<Integer, Integer> baseCounts;
    private GenomeLoc phasingLocus;

    public MultipleBaseCounts(GenomeLoc phasingLoc) {
        this.baseCounts = new HashMap<Integer, Integer>();
        this.phasingLocus = phasingLoc;
    }

    public boolean samePhasingLocAs(GenomeLoc loc) {
        return phasingLocus.equals(loc);
    }

    public void incrementBaseCount(byte base) {
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
        Integer cnt = baseCounts.get(baseIndex);
        if (cnt == null)
            cnt = 0;

        baseCounts.put(baseIndex, cnt + 1);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        sb.append("Base counts");
        for (Map.Entry<Integer, Integer> baseCountEntry : baseCounts.entrySet()) {
            byte base = BaseUtils.baseIndexToSimpleBase(baseCountEntry.getKey());
            int cnt = baseCountEntry.getValue();
            sb.append("\t" + (char) base + ": " + cnt);
        }

        return sb.toString();
    }
}
