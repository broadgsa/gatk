/*
 * Copyright (c) 2010 The Broad Institute
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
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;


import java.util.ArrayList;
import java.util.List;
import java.util.Collection;
import java.io.PrintStream;


/**
 * Given a set of reads and variant data from an external source, this walker downsamples the reads at variant
 * positions to empirically assess the rate at which variants would be confidently and correctly called given different levels of coverage.
 */
public class SnpCallRateByCoverageWalker extends LocusWalker<List<String>, String> {
    @Output
    PrintStream out;

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName="min_confidence_threshold", shortName="confidence", doc="The phred-scaled confidence threshold by which variants should be filtered", required=false) public int confidence = 50;
    @Argument(fullName="min_coverage", shortName="mincov", doc="Mininum coverage to downsample to", required=false) public int min_coverage=1;
    @Argument(fullName="max_coverage", shortName="maxcov", doc="Maximum coverage to downsample to", required=false) public int max_coverage=Integer.MAX_VALUE;
    @Argument(fullName="downsampling_repeats", shortName="repeat", doc="Number of times to repeat downsampling at each coverage level", required=false) public int downsampling_repeats=1;
    @Argument(fullName="coverage_step_size", shortName="step", doc="Coverage step size", required=false) public int step=1;

    UnifiedGenotyperEngine UG;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.STANDARD_CONFIDENCE_FOR_CALLING = uac.STANDARD_CONFIDENCE_FOR_EMITTING = confidence;
        uac.ALL_BASES_MODE = true;
        UG = new UnifiedGenotyperEngine(getToolkit(), uac);

        out.println("#locus\tid\tdownsampled_coverage\tpct_coverage\titeration\tref\teval_call\tcomp_call\tvariant_concordance\tgenotype_concordance");
    }

    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 &&
                context.getBasePileup().size() != 0 &&
                tracker != null &&
                tracker.getAllVariantContexts(ref) != null
        );
    }

    public List<String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Collection<VariantContext> contexts = tracker.getAllVariantContexts(ref);

        for (VariantContext vc : contexts) {
            if (vc.isVariant() && !vc.isFiltered()) {
                //out.println(vc.toString());

                ArrayList<String> GenotypeCalls = new ArrayList<String>();

                List<SAMRecord> reads = context.getReads();
                List<Integer> offsets = context.getOffsets();

                int coverage_available = reads.size();
                List<Integer> coverage_levels = new ArrayList<Integer>();
                //Integer this_max_coverage = Math.min(max_coverage, coverage_available);
                Integer this_max_coverage = 100;
                //for (int coverage = min_coverage; coverage <= this_max_coverage; coverage++) {
                for (int coverage = min_coverage; coverage <= this_max_coverage; coverage += step) {
                    coverage_levels.add(coverage);
                }

                // Iterate over coverage levels
                for (int coverage : coverage_levels) {
                    int usableCoverage = Math.min(coverage_available, coverage); // don't exceed max available coverage

                    Genotype vcCall = vc.getGenotype(0);
                    Genotype call = null;
                    int goodIterations = 0;

                    for (int r=0; r < downsampling_repeats; r++) {
                        List<Integer> subset_indices = MathUtils.sampleIndicesWithReplacement(coverage_available, usableCoverage);
                        List<SAMRecord> sub_reads = MathUtils.sliceListByIndices(subset_indices, reads);
                        List<Integer> sub_offsets = MathUtils.sliceListByIndices(subset_indices, offsets);

                        AlignmentContext subContext = new AlignmentContext(context.getLocation(), new ReadBackedPileupImpl(context.getLocation(),sub_reads, sub_offsets));

                        VariantCallContext calls = UG.runGenotyper(tracker, ref, subContext);

                        if (calls != null && calls.vc != null && calls.vc.getNSamples() > 0 && calls.confidentlyCalled) {
                            Genotype evCall = calls.vc.getGenotype(0);
                            vcCall = vc.getGenotype(evCall.getSampleName());

                            if ((evCall.isHet() || evCall.isHomVar()) && (vcCall.isHet() || vcCall.isHomVar())) {
                                call = evCall;
                                goodIterations++;
                            }

                        }
                    }

                    out.printf("%s\t%s\t\t%d\t%f\t%d\t%c\t%s\t%s\t%d\t%d%n",
                               context.getLocation(),
                               vc.hasAttribute(VariantContext.ID_KEY) ? vc.getAttribute(VariantContext.ID_KEY) : "?",
                               coverage,
                               ((float) coverage)/((float) reads.size()),
                               goodIterations,
                               (char)BaseUtils.baseIndexToSimpleBase(ref.getBaseIndex()),
                               call == null ? "./." : call.getGenotypeString(),
                               vcCall.getGenotypeString(),
                               call == null ? 0 : call.getType() == vcCall.getType() ? 1 : 0,
                               call == null ? 0 : (call.isHet() || call.isHomVar()) && (vcCall.isHet() || vcCall.isHomVar()) ? 1 : 0);
                }
                return GenotypeCalls;
            }
        }
        
        return null;
    }

    public String reduceInit() {
        return "";
    }

    public void onTraversalDone(String result) {} // Don't print the reduce result

    public String reduce(List<String> alleleFreqLines, String sum) {
        /*
        for (String line : alleleFreqLines) {
            out.println(line);
        }

        */
        return "";
	}
}
