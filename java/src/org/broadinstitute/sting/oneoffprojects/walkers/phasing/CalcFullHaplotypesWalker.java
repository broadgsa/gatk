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

package org.broadinstitute.sting.oneoffprojects.walkers.phasing;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.phasing.ReadBackedPhasingWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.PrintStream;
import java.util.*;

import static org.broadinstitute.sting.utils.vcf.VCFUtils.getVCFHeadersFromRods;

/**
 * Walks along all variant ROD loci and verifies the phasing from the reads for user-defined pairs of sites.
 */
@Allows(value = {DataSource.REFERENCE})
@Requires(value = {DataSource.REFERENCE})

public class CalcFullHaplotypesWalker extends RodWalker<Integer, Integer> {
    private Map<String, Haplotype> waitingHaplotypes = null;

    @Output(doc = "File to which results should be written", required = true)
    protected PrintStream out;

    @Argument(doc = "sample to emit", required = false)
    protected String sample = null;

    @Argument(doc = "only include physically-phased results", required = false)
    protected boolean requirePQ = false;

    public void initialize() {
        this.waitingHaplotypes = new HashMap<String, Haplotype>();

        Map<String, VCFHeader> rodNameToHeader = getVCFHeadersFromRods(getToolkit(), null);
        for (VCFHeader header : rodNameToHeader.values()) {
            for (String sample : header.getGenotypeSamples())
                waitingHaplotypes.put(sample, null);
        }
    }

    public boolean generateExtendedEvents() {
        return false;
    }

    public Integer reduceInit() {
        return 0;
    }

    /**
     * @param tracker the meta-data tracker
     * @param ref     the reference base
     * @param context the context for the given locus
     * @return statistics of and list of all phased VariantContexts and their base pileup that have gone out of cacheWindow range.
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null)
            return null;

        GenomeLoc curLocus = ref.getLocus();
        outputDoneHaplotypes(curLocus);

        int curPosition = curLocus.getStop();
        int prevPosition = curPosition - 1;

        // Extend the haplotypes to include up to this position (BUT EXCLUSIVE OF THIS POSITION):
        for (Map.Entry<String, Haplotype> sampleHapEntry : waitingHaplotypes.entrySet()) {
            Haplotype waitingHaplotype = sampleHapEntry.getValue();

            if (waitingHaplotype == null) {// changed to a new contig:
                // Set the new haplotype to extend from [1, prevPosition]
                if (prevPosition >= 1) {
                    GenomeLoc startInterval = getToolkit().getGenomeLocParser().parseGenomeLoc(curLocus.getContig(), 1, prevPosition);
                    waitingHaplotype = new Haplotype(startInterval, sampleHapEntry.getKey());
                    sampleHapEntry.setValue(waitingHaplotype);
                }
            }
            else
                waitingHaplotype.extend(prevPosition);
        }

        Collection<VariantContext> vcs = tracker.getAllVariantContexts(ref, context.getLocation());
        for (VariantContext vc : vcs) {
            if (vc.isFiltered())
                continue;

            if (sample != null)
                vc = vc.subContextFromGenotypes(vc.getGenotype(sample));

            for (Map.Entry<String, Genotype> sampleGtEntry : vc.getGenotypes().entrySet()) {
                String sample = sampleGtEntry.getKey();
                Genotype gt = sampleGtEntry.getValue();

                if (gt.isHet()) {
                    Haplotype sampleHap = waitingHaplotypes.get(sample);
                    if (sampleHap == null)
                        throw new ReviewedStingException("EVERY sample should have a haplotype [by code above and getToolkit().getSamples()]");

                    // Terminate the haplotype before here:
                    if (!gt.isPhased() || (requirePQ && !gt.hasAttribute(ReadBackedPhasingWalker.PQ_KEY))) {
                        outputHaplotype(sampleHap);

                        // Start a new haplotype from the current position:
                        sampleHap = new Haplotype(curLocus, sample);
                        waitingHaplotypes.put(sample, sampleHap);
                    }
                    else {
                        sampleHap.extend(curPosition);
                    }

                    sampleHap.incrementHetCount();
                }
            }
        }

        return 1;
    }

    public Integer reduce(Integer addIn, Integer runningCount) {
        if (addIn == null)
            addIn = 0;

        return runningCount + addIn;
    }

    private void outputDoneHaplotypes(GenomeLoc curLocus) {
        for (Map.Entry<String, Haplotype> sampleHapEntry : waitingHaplotypes.entrySet()) {
            Haplotype waitingHaplotype = sampleHapEntry.getValue();

            if (waitingHaplotype != null) {
                if (curLocus == null || !waitingHaplotype.interval.onSameContig(curLocus)) {
                    sampleHapEntry.setValue(null);

                    // Set the output haplotype to terminate at the end of its contig:
                    int contigLength = getContigLength(waitingHaplotype.interval.getContig());
                    waitingHaplotype.extend(contigLength);
                    outputHaplotype(waitingHaplotype);
                }
            }
        }
    }

    private int getContigLength(String contig) {
        return getToolkit().getGenomeLocParser().getContigInfo(contig).getSequenceLength();
    }

    private void outputHaplotype(Haplotype h) {
        out.println(h);
    }

    /**
     * @param result the number of reads and VariantContexts seen.
     */
    public void onTraversalDone(Integer result) {
        outputDoneHaplotypes(null);

        System.out.println("map was called " + result + " times.");
    }

    private class Haplotype {
        public GenomeLoc interval;
        public String sample;
        public int hetCount;

        public Haplotype(GenomeLoc interval, String sample) {
            this.interval = interval;
            this.sample = sample;
            this.hetCount = 0;
        }

        public void extend(int stop) {
            if (stop > interval.getStop())
                interval = getToolkit().getGenomeLocParser().parseGenomeLoc(interval.getContig(), interval.getStart(), stop);
        }

        public void incrementHetCount() {
            hetCount++;
        }

        public String toString() {
            return sample + "\t" + interval.toString() + "\t" + interval.size() + "\t" + hetCount;
        }
    }
}