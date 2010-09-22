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

package org.broadinstitute.sting.playground.gatk.walkers.phasing;

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
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import javax.naming.OperationNotSupportedException;
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
    private DoublyLinkedList<UnfinishedVariantAndReads> partiallyPhasedSites = null; // the phased VCs to be emitted, and the alignment bases at these positions

    private static PreciseNonNegativeDouble ZERO = new PreciseNonNegativeDouble(0.0);

    private LinkedList<String> rodNames = null;
    private PhasingQualityStatsWriter statsWriter = null;

    public void initialize() {
        if (maxPhaseSites <= 2)
            maxPhaseSites = 2; // by definition, must phase a site relative to previous site [thus, 2 in total]

        rodNames = new LinkedList<String>();
        rodNames.add("variant");

        unphasedSiteQueue = new LinkedList<VariantAndReads>();
        partiallyPhasedSites = new DoublyLinkedList<UnfinishedVariantAndReads>();

        initializeVcfWriter();

        if (variantStatsFilePrefix != null)
            statsWriter = new PhasingQualityStatsWriter(variantStatsFilePrefix);
    }

    private void initializeVcfWriter() {
        // setup the header fields:
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
                          (note that we ASSUME that the VCF is ordered by <contig,locus>).
                           Note that this will always leave at least one entry (the last one), since lastLocus is in range of itself.
                         */
                        break;
                    }
                    // Already saw all variant positions within cacheWindow distance ahead of vc (on its contig)
                }
                // Update partiallyPhasedSites before it's used in phaseSite:
                oldPhasedList.addAll(discardIrrelevantPhasedSites());
                logger.debug("oldPhasedList(1st) = " + toStringVCL(oldPhasedList));

                VariantAndReads vr = unphasedSiteQueue.remove();
                logger.debug("Performing phasing for " + VariantContextUtils.getLocation(vr.variant));
                phaseSite(vr, phaseStats);
            }
        }

        // Update partiallyPhasedSites after phaseSite is done:
        oldPhasedList.addAll(discardIrrelevantPhasedSites());
        logger.debug("oldPhasedList(2nd) = " + toStringVCL(oldPhasedList));
        return oldPhasedList;
    }

    private List<VariantContext> discardIrrelevantPhasedSites() {
        List<VariantContext> vcList = new LinkedList<VariantContext>();

        GenomeLoc nextToPhaseLoc = null;
        if (!unphasedSiteQueue.isEmpty())
            nextToPhaseLoc = VariantContextUtils.getLocation(unphasedSiteQueue.peek().variant);

        while (!partiallyPhasedSites.isEmpty()) {
            if (nextToPhaseLoc != null) { // otherwise, unphasedSiteQueue.isEmpty(), and therefore no need to keep any of the "past"
                UnfinishedVariantAndReads partPhasedVr = partiallyPhasedSites.peek();

                if (partPhasedVr.processVariant && isInWindowRange(partPhasedVr.unfinishedVariant.getLocation(), nextToPhaseLoc))
                    // nextToPhaseLoc is still not far enough ahead of partPhasedVr to exclude partPhasedVr from calculations
                    break;
            }
            vcList.add(partiallyPhasedSites.remove().unfinishedVariant.toVariantContext());
        }

        return vcList;
    }

    /* Phase vc (removed head of unphasedSiteQueue) using all VariantContext objects in
       partiallyPhasedSites, and all in unphasedSiteQueue that are within cacheWindow distance ahead of vc (on its contig).

       ASSUMES: All VariantContexts in unphasedSiteQueue are in positions downstream of vc (head of queue).
     */
    private void phaseSite(VariantAndReads vr, PhasingStats phaseStats) {
        UnfinishedVariantAndReads pvr = new UnfinishedVariantAndReads(vr);
        if (!vr.processVariant) {
            partiallyPhasedSites.add(pvr);
            return;
        }

        VariantContext vc = vr.variant;
        logger.debug("Will phase vc = " + VariantContextUtils.getLocation(vc));
        UnfinishedVariantContext uvc = pvr.unfinishedVariant;

        // Perform per-sample phasing:
        Map<String, Genotype> sampGenotypes = vc.getGenotypes();
        Map<String, PhaseCounts> samplePhaseStats = new TreeMap<String, PhaseCounts>();
        for (Map.Entry<String, Genotype> sampGtEntry : sampGenotypes.entrySet()) {
            String samp = sampGtEntry.getKey();
            Genotype gt = sampGtEntry.getValue();

            logger.debug("sample = " + samp);
            if (isCalledDiploidGenotype(gt) && gt.isHet()) { // Can attempt to phase this genotype
                PhasingWindow phaseWindow = new PhasingWindow(vr, samp);
                if (phaseWindow.hasPreviousHets()) { // Otherwise, nothing to phase this against
                    BialleleSNP biall = new BialleleSNP(gt);
                    logger.debug("Want to phase TOP vs. BOTTOM for: " + "\n" + biall);

                    DoublyLinkedList.BidirectionalIterator<UnfinishedVariantAndReads> prevHetAndInteriorIt = phaseWindow.prevHetAndInteriorIt;
                    /* Notes:
                     1. Call to next() advances iterator to next position in partiallyPhasedSites.
                     2. prevHetGenotype != null, since otherwise prevHetAndInteriorIt would not have been chosen to point to its UnfinishedVariantAndReads.
                     */
                    UnfinishedVariantContext prevUvc = prevHetAndInteriorIt.next().unfinishedVariant;
                    Genotype prevHetGenotype = prevUvc.getGenotype(samp);

                    PhaseResult pr = phaseSample(phaseWindow);
                    boolean genotypesArePhased = (pr.phaseQuality >= phaseQualityThresh);
                    if (genotypesArePhased) {
                        BialleleSNP prevBiall = new BialleleSNP(prevHetGenotype);

                        logger.debug("THE PHASE PREVIOUSLY CHOSEN FOR PREVIOUS:\n" + prevBiall + "\n");
                        logger.debug("THE PHASE CHOSEN HERE:\n" + biall + "\n\n");

                        ensurePhasing(biall, prevBiall, pr.haplotype);
                        Map<String, Object> gtAttribs = new HashMap<String, Object>(gt.getAttributes());
                        gtAttribs.put("PQ", pr.phaseQuality);
                        Genotype phasedGt = new Genotype(gt.getSampleName(), biall.getAllelesAsList(), gt.getNegLog10PError(), gt.getFilters(), gtAttribs, genotypesArePhased);
                        uvc.setGenotype(samp, phasedGt);
                    }

                    // Now, update the 0 or more "interior" hom sites in between the previous het site and this het site:
                    while (prevHetAndInteriorIt.hasNext()) {
                        UnfinishedVariantAndReads interiorVr = prevHetAndInteriorIt.next();
                        if (interiorVr.processVariant) {
                            UnfinishedVariantContext interiorUvc = interiorVr.unfinishedVariant;
                            Genotype handledGt = interiorUvc.getGenotype(samp);
                            if (handledGt == null || !isCalledDiploidGenotype(handledGt))
                                throw new ReviewedStingException("LOGICAL error: should not have breaks WITHIN haplotype");
                            if (!handledGt.isHom())
                                throw new ReviewedStingException("LOGICAL error: should not have anything besides hom sites IN BETWEEN two het sites");

                            // Use the same PQ for each hom site in the "interior" as for the het-het phase:
                            if (genotypesArePhased) {
                                Map<String, Object> handledGtAttribs = new HashMap<String, Object>(handledGt.getAttributes());
                                handledGtAttribs.put("PQ", pr.phaseQuality);
                                Genotype phasedHomGt = new Genotype(handledGt.getSampleName(), handledGt.getAlleles(), handledGt.getNegLog10PError(), handledGt.getFilters(), handledGtAttribs, genotypesArePhased);
                                interiorUvc.setGenotype(samp, phasedHomGt);
                            }
                        }
                    }

                    if (statsWriter != null)
                        statsWriter.addStat(samp, VariantContextUtils.getLocation(vc), distance(prevUvc, vc), pr.phaseQuality, phaseWindow.readsAtHetSites.size(), phaseWindow.hetGenotypes.length);

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
        }

        partiallyPhasedSites.add(pvr); // only add it in now, since don't want it to be there during phasing
        phaseStats.addIn(new PhasingStats(samplePhaseStats));
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
        private DoublyLinkedList.BidirectionalIterator<UnfinishedVariantAndReads> prevHetAndInteriorIt = null;
        private int phasingSiteIndex = -1;
        private Map<String, Read> readsAtHetSites = null;

        public boolean hasPreviousHets() {
            return phasingSiteIndex > 0;
        }

        // ASSUMES that: isCalledDiploidGenotype(gt) && gt.isHet() [gt = vr.unfinishedVariant.getGenotype(sample)]
        public PhasingWindow(VariantAndReads vr, String sample) {
            List<GenotypeAndReadBases> listHetGenotypes = new LinkedList<GenotypeAndReadBases>();

            // Include previously phased sites in the phasing computation:
            DoublyLinkedList.BidirectionalIterator<UnfinishedVariantAndReads> phasedIt = partiallyPhasedSites.iterator();
            while (phasedIt.hasNext()) {
                UnfinishedVariantAndReads phasedVr = phasedIt.next();
                if (phasedVr.processVariant) {
                    Genotype gt = phasedVr.unfinishedVariant.getGenotype(sample);
                    if (gt == null || !isCalledDiploidGenotype(gt)) { // constructed haplotype must start AFTER this "break"
                        listHetGenotypes.clear(); // clear out any history
                    }
                    else if (gt.isHet()) {
                        GenotypeAndReadBases grb = new GenotypeAndReadBases(gt, phasedVr.sampleReadBases.get(sample), phasedVr.unfinishedVariant.getLocation());
                        listHetGenotypes.add(grb);
                        logger.debug("Using UPSTREAM het site = " + grb.loc);
                        prevHetAndInteriorIt = phasedIt.clone();
                    }
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
            GenomeLoc phaseLocus = VariantContextUtils.getLocation(vr.variant);
            GenotypeAndReadBases grbPhase = new GenotypeAndReadBases(vr.variant.getGenotype(sample), vr.sampleReadBases.get(sample), phaseLocus);
            listHetGenotypes.add(grbPhase);
            logger.debug("PHASING het site = " + grbPhase.loc + " [phasingSiteIndex = " + phasingSiteIndex + "]");

            // Include as-of-yet unphased sites in the phasing computation:
            for (VariantAndReads nextVr : unphasedSiteQueue) {
                if (!isInWindowRange(vr.variant, nextVr.variant)) //nextVr too far ahead of the range used for phasing vc
                    break;
                if (nextVr.processVariant) {
                    Genotype gt = nextVr.variant.getGenotype(sample);
                    if (gt == null || !isCalledDiploidGenotype(gt)) { // constructed haplotype must end BEFORE this "break"
                        break;
                    }
                    else if (gt.isHet()) {
                        GenotypeAndReadBases grb = new GenotypeAndReadBases(gt, nextVr.sampleReadBases.get(sample), VariantContextUtils.getLocation(nextVr.variant));
                        listHetGenotypes.add(grb);
                        logger.debug("Using DOWNSTREAM het site = " + grb.loc);
                    }
                }
            }

            // First, assemble the "sub-reads" from the COMPLETE WINDOW-BASED SET of heterozygous positions for this sample:
            buildReadsAtHetSites(listHetGenotypes);

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
            hetGenotypes = new Genotype[listHetGenotypes.size()];
            int index = 0;
            for (GenotypeAndReadBases copyGrb : listHetGenotypes)
                hetGenotypes[index++] = copyGrb.genotype;
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes) {
            buildReadsAtHetSites(listHetGenotypes, null);
        }

        private void buildReadsAtHetSites(List<GenotypeAndReadBases> listHetGenotypes, Set<String> onlyKeepReads) {
            readsAtHetSites = new HashMap<String, Read>();

            LinkedList<ReadBasesAtPosition> basesAtPositions = new LinkedList<ReadBasesAtPosition>();
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                ReadBasesAtPosition readBases = grb.readBases;
                if (readBases == null)
                    readBases = new ReadBasesAtPosition(); // for transparency, put an empty list of bases at this position for sample
                basesAtPositions.add(readBases);
            }

            int index = 0;
            for (ReadBasesAtPosition rbp : basesAtPositions) {
                for (ReadBase rb : rbp) {
                    String readName = rb.readName;
                    if (onlyKeepReads != null && !onlyKeepReads.contains(readName)) // if onlyKeepReads exists, ignore reads not in onlyKeepReads
                        continue;

                    Read rd = readsAtHetSites.get(readName);
                    if (rd == null) {
                        rd = new Read(basesAtPositions.size(), rb.mappingQual);
                        readsAtHetSites.put(readName, rd);
                    }
                    rd.updateBaseAndQuality(index, rb.base, rb.baseQual);
                }
                index++;
            }
            logger.debug("Number of sites in window = " + index);

            if (logger.isDebugEnabled()) {
                logger.debug("ALL READS:");
                for (Map.Entry<String, Read> nameToReads : readsAtHetSites.entrySet()) {
                    String rdName = nameToReads.getKey();
                    Read rd = nameToReads.getValue();
                    logger.debug(rd + "\t" + rdName);
                }
            }
        }

        private class ReadProperties {
            public List<GraphEdge> rdEdges;
            public int[] siteInds;

            public ReadProperties(Read rd) {
                this.siteInds = rd.getNonNullIndices();
                this.rdEdges = new LinkedList<GraphEdge>();

                // sufficient to create a path linking the sites in rd, so they all end up in the same connected component:
                for (int i = 0; i < siteInds.length - 1; i++) {
                    GraphEdge e = new GraphEdge(siteInds[i], siteInds[i + 1]);
                    rdEdges.add(e);
                }
            }
        }

        private class EdgeCounts {
            private Map<GraphEdge, Integer> counts;

            public EdgeCounts() {
                this.counts = new TreeMap<GraphEdge, Integer>(); // implemented GraphEdge.compareTo()
            }

            public int getCount(GraphEdge e) {
                Integer count = counts.get(e);
                if (count == null)
                    return 0;

                return count;
            }

            public int incrementEdge(GraphEdge e) {
                Integer eCount = counts.get(e);
                int cnt;
                if (eCount == null)
                    cnt = 0;
                else
                    cnt = eCount;

                cnt++;
                counts.put(e, cnt);
                return cnt;
            }

            public int decrementEdge(GraphEdge e) {
                Integer eCount = counts.get(e);
                if (eCount == null)
                    return 0;

                int cnt = eCount - 1;
                counts.put(e, cnt);
                return cnt;
            }
        }

        public Set<String> removeExtraneousReads(int numHetSites) {
            Graph readGraph = new Graph(numHetSites);
            Map<String, ReadProperties> readToGraphProperties = new HashMap<String, ReadProperties>();
            EdgeCounts edgeCounts = new EdgeCounts();

            for (Map.Entry<String, Read> nameToReads : readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                Read rd = nameToReads.getValue();

                ReadProperties rp = new ReadProperties(rd);
                if (!rp.rdEdges.isEmpty()) { // otherwise, this read is clearly irrelevant since it can't link anything
                    for (GraphEdge e : rp.rdEdges) {
                        readGraph.addEdge(e);
                        logger.debug("Read = " + rdName + " is adding edge: " + e);

                        edgeCounts.incrementEdge(e);
                    }
                    readToGraphProperties.put(rdName, rp);
                }
            }
            logger.debug("Read graph:\n" + readGraph);
            Set<String> keepReads = new HashSet<String>();

            // Check which Reads are involved in paths from (phasingSiteIndex - 1) to (phasingSiteIndex):
            int prev = phasingSiteIndex - 1;
            int cur = phasingSiteIndex;

            if (!readGraph.getConnectedComponents().inSameSet(prev, cur)) { // There is NO path between cur and prev
                logger.debug("NO READ PATH between PHASE site [" + cur + "] and UPSTREAM site [" + prev + "]");
                readsAtHetSites.clear();
                return keepReads;
            }

            for (Map.Entry<String, ReadProperties> rdEdgesEntry : readToGraphProperties.entrySet()) {
                String testRead = rdEdgesEntry.getKey();
                ReadProperties rp = rdEdgesEntry.getValue();
                logger.debug("Testing the connectivity of Read: " + testRead);

                // Check the connected components after removing this read's UNIQUE edges:
                for (GraphEdge e : rp.rdEdges) {
                    if (edgeCounts.getCount(e) == 1) // otherwise, the edge still exists without this read
                        readGraph.removeEdge(e);
                }
                DisjointSet ccAfterRemove = readGraph.getConnectedComponents();

                /* testRead contributes a path between prev and cur iff:
                   There exists i != j s.t. testRead[i] != null, testRead[j] != null, ccAfterRemove.inSameSet(prev,i) && ccAfterRemove.inSameSet(j,cur)
                   [since ALL non-null indices in testRead are connected to one another, as one clique].
                 */
                List<Integer> sameCCasPrev = ccAfterRemove.inSameSetAs(prev, rp.siteInds);
                List<Integer> sameCCasCur = ccAfterRemove.inSameSetAs(cur, rp.siteInds);
                if (logger.isDebugEnabled()) {
                    StringBuilder sb = new StringBuilder("sameCCasPrev:");
                    for (int ind : sameCCasPrev)
                        sb.append(" " + ind);
                    logger.debug(sb.toString());

                    sb = new StringBuilder("sameCCasCur:");
                    for (int ind : sameCCasCur)
                        sb.append(" " + ind);
                    logger.debug(sb.toString());
                }

                boolean keepRead = false;
                if (!sameCCasPrev.isEmpty() && !sameCCasCur.isEmpty()) { // There exists a path from prev to cur that goes through the sites in testRead
                    // Now, make sure that TWO DISTINCT sites, i and j, in testRead are used in the path:
                    Set<Integer> union = new HashSet<Integer>(sameCCasPrev);
                    union.addAll(sameCCasCur);
                    if (union.size() >= 2) // i != j
                        keepRead = true;
                }

                if (keepRead) {
                    logger.debug("Read is part of path from " + prev + " to " + cur);
                    keepReads.add(testRead);

                    // Add the removed edges back in, since we're keeping the read:
                    for (GraphEdge e : rp.rdEdges)
                        readGraph.addEdge(e);
                }
                else { // Decrease the count for the edges [note that any read-specific edges were already removed above]:
                    for (GraphEdge e : rp.rdEdges)
                        edgeCounts.decrementEdge(e);
                }
            }

            // Retain only the reads that contain an edge in a path connecting prev and cur:
            Iterator<Map.Entry<String, Read>> readIt = readsAtHetSites.entrySet().iterator();
            while (readIt.hasNext()) {
                Map.Entry<String, Read> nameToReads = readIt.next();
                String rdName = nameToReads.getKey();
                if (!keepReads.contains(rdName)) {
                    readIt.remove();
                    logger.debug("Removing extraneous read: " + rdName);
                }
            }

            return keepReads;
        }
        
        private List<GenotypeAndReadBases> removeExtraneousSites(List<GenotypeAndReadBases> listHetGenotypes) {
            Set<Integer> sitesWithReads = new HashSet<Integer>();
            for (Map.Entry<String, Read> nameToReads : readsAtHetSites.entrySet()) {
                Read rd = nameToReads.getValue();
                for (int i : rd.getNonNullIndices())
                    sitesWithReads.add(i);
            }

            // Remove all sites that have no read bases:
            List<GenotypeAndReadBases> keepHetSites = new LinkedList<GenotypeAndReadBases>();
            int index = 0;
            int numPrecedingRemoved = 0;
            for (GenotypeAndReadBases grb : listHetGenotypes) {
                boolean keepSite = sitesWithReads.contains(index);
                if (logger.isDebugEnabled() && !keepSite)
                    logger.debug("Removing read-less site " + grb.loc);

                if (keepSite || index == phasingSiteIndex || index == phasingSiteIndex - 1) {
                    keepHetSites.add(grb);
                    if (!keepSite)
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
            logger.warn("REDUCED to " + listHetGenotypes.size() + " sites:\n" + toStringGRL(listHetGenotypes));

            return listHetGenotypes;
        }
    }

    private PhaseResult phaseSample(PhasingWindow phaseWindow) {
        /* Will map a phase and its "complement" to a single representative phase,
          and marginalizeTable() marginalizes to 2 positions [starting at the previous position, and then the current position]:
        */
        HaplotypeTableCreator tabCreator = new BiallelicComplementHaplotypeTableCreator(phaseWindow.hetGenotypes, phaseWindow.phasingSiteIndex - 1, 2);
        PhasingTable sampleHaps = tabCreator.getNewTable();

        logger.debug("Number of USED reads [connecting the two positions to be phased] at sites: " + phaseWindow.readsAtHetSites.size());
        if (logger.isDebugEnabled()) {
            logger.debug("USED READS:");
            for (Map.Entry<String, Read> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
                String rdName = nameToReads.getKey();
                Read rd = nameToReads.getValue();
                logger.debug(rd + "\t" + rdName);
            }
        }

        // Update the phasing table based on each of the sub-reads for this sample:
        for (Map.Entry<String, Read> nameToReads : phaseWindow.readsAtHetSites.entrySet()) {
            Read rd = nameToReads.getValue();

            logger.debug("rd = " + rd + "\tname = " + nameToReads.getKey() + (rd.isGapped() ? "\tGAPPED" : ""));

            for (PhasingTable.PhasingTableEntry pte : sampleHaps) {
                PhasingScore score = rd.matchHaplotypeClassScore(pte.getHaplotypeClass());
                pte.getScore().integrateReadScore(score);

                logger.debug("score(" + rd + ", " + pte.getHaplotypeClass() + ") = " + score);
            }
        }
        logger.debug("\nPhasing table [AFTER CALCULATION]:\n" + sampleHaps + "\n");

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

        return new PhaseResult(maxEntry.getHaplotypeClass().getRepresentative(), phaseQuality);
    }

    /*
        Ensure that curBiall is phased relative to prevBiall as specified by hap.
     */

    public static void ensurePhasing(BialleleSNP curBiall, BialleleSNP prevBiall, Haplotype hap) {
        if (hap.size() < 2)
            throw new ReviewedStingException("LOGICAL ERROR: Only considering haplotypes of length > 2!");

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

    private static int distance(GenomeLoc loc1, GenomeLoc loc2) {
        if (!loc1.onSameContig(loc2))
            return Integer.MAX_VALUE;

        return loc1.distance(loc2);
    }

    private static int distance(VariantContext vc1, VariantContext vc2) {
        GenomeLoc loc1 = VariantContextUtils.getLocation(vc1);
        GenomeLoc loc2 = VariantContextUtils.getLocation(vc2);

        return distance(loc1, loc2);
    }

    private static int distance(UnfinishedVariantContext uvc1, VariantContext vc2) {
        GenomeLoc loc1 = uvc1.getLocation();
        GenomeLoc loc2 = VariantContextUtils.getLocation(vc2);

        return distance(loc1, loc2);
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

        System.out.println("\n-- Phasing summary [minimal haplotype quality (PQ): " + phaseQualityThresh + "] --");
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

    private static class UnfinishedVariantAndReads {
        public UnfinishedVariantContext unfinishedVariant;
        public HashMap<String, ReadBasesAtPosition> sampleReadBases;
        public boolean processVariant;

        public UnfinishedVariantAndReads(UnfinishedVariantContext unfinishedVariant, HashMap<String, ReadBasesAtPosition> sampleReadBases, boolean processVariant) {
            this.unfinishedVariant = unfinishedVariant;
            this.sampleReadBases = sampleReadBases;
            this.processVariant = processVariant;
        }

        public UnfinishedVariantAndReads(VariantAndReads vr) {
            this.unfinishedVariant = new UnfinishedVariantContext(vr.variant);
            this.sampleReadBases = vr.sampleReadBases;
            this.processVariant = vr.processVariant;
        }
    }

    // COULD replace with MutableVariantContext if it worked [didn't throw exceptions when trying to call its set() methods]...
    private static class UnfinishedVariantContext {
        private String name;
        private String contig;
        private long start;
        private long stop;
        private Collection<Allele> alleles;
        private Map<String, Genotype> genotypes;
        private double negLog10PError;
        private Set<String> filters;
        private Map<String, ?> attributes;

        public UnfinishedVariantContext(VariantContext vc) {
            this.name = vc.getName();
            this.contig = vc.getChr();
            this.start = vc.getStart();
            this.stop = vc.getEnd();
            this.alleles = vc.getAlleles();
            this.genotypes = new HashMap<String, Genotype>(vc.getGenotypes()); // since vc.getGenotypes() is unmodifiable
            this.negLog10PError = vc.getNegLog10PError();
            this.filters = vc.getFilters();
            this.attributes = vc.getAttributes();
        }

        public VariantContext toVariantContext() {
            return new VariantContext(name, contig, start, stop, alleles, genotypes, negLog10PError, filters, attributes);
        }

        public GenomeLoc getLocation() {
            return GenomeLocParser.createGenomeLoc(contig, start, stop);
        }

        public Genotype getGenotype(String sample) {
            return genotypes.get(sample);
        }

        public void setGenotype(String sample, Genotype newGt) {
            genotypes.put(sample, newGt);
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

        public boolean isEmpty() {
            return bases.isEmpty();
        }
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
                    hapBases[i] = BialleleSNP.getSingleBase(alleleI);
                }
                allHaps.add(new Haplotype(hapBases));
            }
            return allHaps;
        }

        public static PhasingTable marginalizeTable(PhasingTable table) {
            Map<Haplotype, PreciseNonNegativeDouble> hapMap = new TreeMap<Haplotype, PreciseNonNegativeDouble>();
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

        public BiallelicComplementHaplotypeTableCreator(Genotype[] hetGenotypes, int startIndex, int marginalizeLength) {
            super(hetGenotypes);

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
                throw new ReviewedStingException("INTERNAL ERROR: hap.size() != numSites");

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

        public PhaseResult(Haplotype haplotype, double phaseQuality) {
            this.haplotype = haplotype;
            this.phaseQuality = phaseQuality;
        }
    }

    public static boolean isUnfilteredBiallelicSNP(VariantContext vc) {
        return (vc.isSNP() && vc.isBiallelic() && !vc.isFiltered());
    }

    public static boolean isCalledDiploidGenotype(Genotype gt) {
        return (gt.isCalled() && gt.getPloidy() == 2);
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
            throw new ReviewedStingException("Internal error: CANNOT have null for a missing Haplotype base!");
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
            throw new ReviewedStingException("Read and Haplotype should have same length to be compared!");

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

    public int[] getNonNullIndices() {
        List<Integer> nonNull = new LinkedList<Integer>();
        for (int i = 0; i < bases.length; i++) {
            if (getBase(i) != null)
                nonNull.add(i);
        }

        int[] nonNullArray = new int[nonNull.size()];
        int index = 0;
        for (int i : nonNull)
            nonNullArray[index++] = i;
        return nonNullArray;
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

    public void addStat(String sample, GenomeLoc locus, int distanceFromPrevious, double phasingQuality, int numReads, int windowSize) {
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
            sampWriter.write(locus + "\t" + distanceFromPrevious + "\t" + phasingQuality + "\t" + numReads + "\t" + windowSize + "\n");
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