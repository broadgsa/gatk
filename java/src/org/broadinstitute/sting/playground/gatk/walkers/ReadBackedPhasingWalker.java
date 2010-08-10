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
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.*;
import java.util.*;


/**
 * Walks along all loci, caching a user-defined window of VariantContext sites, and then finishes phasing them when they go out of range (using downstream reads).
 * Use '-BTI variant' to only stop at positions in the VCF file bound to 'variant'.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= ReferenceOrderedDatum.class))
public class ReadBackedPhasingWalker extends LocusWalker<Pair<VariantContextStats, List<VariantContext>>, VariantContextStats> {

    @Argument(fullName="cacheWindowSize", shortName="cacheWindow", doc="The window size (in bases) to cache variant sites and their reads; [default:300]", required=false)
    protected Integer cacheWindow = 300;

    @Argument(fullName="phasedVCFFile", shortName="phasedVCF", doc="The name of the phased VCF file output", required=true)
    protected String phasedVCFFile = null;

    private VCFWriter writer = null;

    private LinkedList<VariantAndAlignment> siteQueue = null;

    private void initializeVcfWriter(VariantContext vc) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        writer = new VCFWriter(new File(phasedVCFFile));
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(vc.getSampleNames())));
    }

    public void initialize() {
	    siteQueue = new LinkedList<VariantAndAlignment>();
    }

    public boolean generateExtendedEvents() { // want to see indels
        return true;
    }

    public VariantContextStats reduceInit() { return new VariantContextStats(); }

    /**
     * For each site of interest, cache the current site and then use the cache to phase all upstream sites
     * for which "sufficient" information has already been observed.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Pair<VariantContextStats, List<VariantContext>> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContextStats vcStats = new VariantContextStats();
        List<VariantContext> phasedList = new LinkedList<VariantContext>();
        if ( tracker == null )
            return new Pair<VariantContextStats, List<VariantContext>>(vcStats, phasedList);

        List<Object> rods = tracker.getReferenceMetaData("variant");
        ListIterator<Object> rodIt = rods.listIterator();
        while (rodIt.hasNext()) {
            VariantContext vc = VariantContextAdaptors.toVariantContext("variant", rodIt.next(), ref);
            if (vc.getType() == VariantContext.Type.MNP) {
                throw new StingException("Doesn't support phasing for multiple-nucleotide polymorphism!");
            }
            VariantAndAlignment va = new VariantAndAlignment(vc, context);
	        siteQueue.add(va);

            int numReads = 0;
            if (context.hasBasePileup()) {
                numReads = context.getBasePileup().size();
            }
            else if (context.hasExtendedEventPileup()) {
                numReads = context.getExtendedEventPileup().size();
            }
            VariantContextStats addInVcStats = new VariantContextStats(numReads, 1);
            vcStats.addTo(addInVcStats);
	    }

        GenomeLoc refLoc = ref.getLocus();
    	while (!siteQueue.isEmpty()) {
    	    VariantContext vc = siteQueue.peek().variant;
            if (!isInWindowRange(refLoc, VariantContextUtils.getLocation(vc))) { // Already saw all variant positions within cacheWindow distance ahead of vc (on its contig)
                VariantContext phasedVc = this.phaseVariantAndRemove();
    		    phasedList.add(phasedVc);
    	    }
    	    else { // refLoc is still not far enough ahead of vc
    		    break; // since we ASSUME that the VCF is ordered by <contig,locus>
    	    }
	    }

        return new Pair<VariantContextStats, List<VariantContext>>(vcStats, phasedList);
    }

    /* Phase vc (head of siteQueue) using all VariantContext objects in the siteQueue that are
        within cacheWindow distance ahead of vc (on its contig).
        ASSUMES:
        1. siteQueue is NOT empty.
        2. All VariantContexts in siteQueue are in positions downstream of vc (head of queue).
     */
    private VariantContext phaseVariantAndRemove() {
        VariantContext vc = siteQueue.peek().variant;

        ListIterator<VariantAndAlignment> windowIt = siteQueue.listIterator();
        int toIndex = 0;
        while (windowIt.hasNext()) {
            if (isInWindowRange(vc, windowIt.next().variant)) {
                toIndex++;
            }
            else { //moved past the relevant range used for phasing
                break;
            }
        }
        List<VariantAndAlignment> windowVcList = siteQueue.subList(0,toIndex);

        //
        if (true) {
            out.println("Will phase vc = " + VariantContextUtils.getLocation(vc));
            ListIterator<VariantAndAlignment> windowVcIt = windowVcList.listIterator();
            while (windowVcIt.hasNext()) {
                VariantContext phaseInfoVc = windowVcIt.next().variant;
                out.println("Using phaseInfoVc = " + VariantContextUtils.getLocation(phaseInfoVc));
            }
            out.println("");
        }
        //

        Map<String, Genotype> sampGenotypes = vc.getGenotypes();
        Map<String, Genotype> phasedGtMap = new TreeMap<String, Genotype>();

        for (Map.Entry<String, Genotype> entry : sampGenotypes.entrySet()) {
            String samp = entry.getKey();
            Genotype gt = entry.getValue();

            if (gt.getPloidy() != 2) {
                throw new StingException("Doesn't support phasing for ploidy that is not 2!");
            }
            Allele topAll = gt.getAllele(0);
            Allele botAll = gt.getAllele(1);

            ListIterator<VariantAndAlignment> windowVcIt = windowVcList.listIterator();
            while (windowVcIt.hasNext()) {
                VariantAndAlignment va = windowVcIt.next();
                VariantContext phaseInfoVc = va.variant;
                AlignmentContext phaseInfoContext = va.alignment;

                ReadBackedPileup reads = null;
                if (phaseInfoContext.hasBasePileup()) {
                    reads = phaseInfoContext.getBasePileup();
                }
                else if (phaseInfoContext.hasExtendedEventPileup()) {
                    reads = phaseInfoContext.getExtendedEventPileup();
                }
                if (reads != null) {
                    ReadBackedPileup sampleReads = null;
                    if (reads.getSamples().contains(samp)) {
                        // Update the phasing table based on the reads for this sample:
                        sampleReads = reads.getPileupForSample(samp);
                        for (PileupElement p : sampleReads) {
                            SAMRecord rd = p.getRead();
                            out.println("read = " + rd);
                        }
                    }
                }
            }

            Random rn = new Random();
            boolean genotypesArePhased = (rn.nextDouble() > 0.5);

            boolean swapChromosomes = (rn.nextDouble() > 0.5);
            if (swapChromosomes) {
                Allele tmp = topAll;
                topAll = botAll;
                botAll = tmp;
            }
            List<Allele> phasedAll = new ArrayList<Allele>();
            phasedAll.add(0, topAll);
            phasedAll.add(1, botAll);

            Genotype phasedGt = new Genotype(gt.getSampleName(), phasedAll, gt.getNegLog10PError(), gt.getFilters(), gt.getAttributes(), genotypesArePhased);
            phasedGtMap.put(samp, phasedGt);
        }
        siteQueue.remove(); // remove vc from head of queue

        return new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), phasedGtMap, vc.getNegLog10PError(), vc.getFilters(), vc.getAttributes());
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
        if ( writer == null )
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

    public VariantContextStats reduce(Pair<VariantContextStats, List<VariantContext>> statsAndList, VariantContextStats stats) {
	    Iterator<VariantContext> varContIter = statsAndList.second.iterator();
	    writeVarContIter(varContIter);

        stats.addTo(statsAndList.first);
        return stats;
    }

    /**
     * Phase anything left in the cached siteQueue, and report the number of reads and VariantContexts processed.
     *
     * @param result  the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(VariantContextStats result) {
        List<VariantContext> finalList = new LinkedList<VariantContext>();
        while (!siteQueue.isEmpty()) {
            VariantContext phasedVc = this.phaseVariantAndRemove();
            finalList.add(phasedVc);
        }
        writeVarContIter(finalList.iterator());

        if ( writer != null )
            writer.close();

        out.println("Number of reads observed: " + result.getNumReads());
	    out.println("Number of variant sites observed: " + result.getNumVarSites());
        out.println("Average coverage: " + ((double) result.getNumReads() / result.getNumVarSites()));
    }

    protected void writeVarContIter(Iterator<VariantContext> varContIter) {
	    while (varContIter.hasNext()) {
	        VariantContext vc = varContIter.next();
	        writeVCF(vc);
	    }
    }

    private static class VariantAndAlignment {
        public VariantContext variant;
        public AlignmentContext alignment;

        public VariantAndAlignment(VariantContext variant, AlignmentContext alignment) {
            this.variant = variant;
            this.alignment = alignment;
        }
    }
}


class VariantContextStats {
    private int numReads;
    private int numVarSites;

    public VariantContextStats() {
        this.numReads = 0;
        this.numVarSites = 0;
    }

    public VariantContextStats(int numReads, int numVarSites) {
        this.numReads = numReads;
        this.numVarSites = numVarSites;
    }

    public void addTo(VariantContextStats other) {
        this.numReads += other.numReads;
        this.numVarSites += other.numVarSites;
    }

    public int getNumReads() {return numReads;}
    public int getNumVarSites() {return numVarSites;}
}