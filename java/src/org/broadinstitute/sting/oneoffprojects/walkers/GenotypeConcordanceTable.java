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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits a table of chrom/pos/sample/GQ/LoDGivenXX for AA, AB, BB/called genotypes for eval and comp for eval and comp VCFs
 */
@Requires(value={})
public class GenotypeConcordanceTable extends RodWalker<Integer, Integer> {
    public static final String EVAL_NAME = "eval";
    public static final String COMP_NAME = "comp";

    @Output(doc="File to which results should be written",required=true)
    protected PrintStream out;

    @Argument(doc="If provided, we will include information where EVAL is missing and COMP is missing", required=false)
    protected boolean keepDoubleMissing = false;

    @Argument(doc="If provided, we will include information where EVAL is no-called", required=false)
    protected boolean keepEvalNoCall = false;

    @Argument(doc="If provided, we will include information where COMP is no-called", required=false)
    protected boolean keepCompNoCall = false;

    Set<String> samples;

    @Override
    public void initialize() {
        samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), Arrays.asList(EVAL_NAME));
        logger.info("Samples: " + samples);
        List<String> fields = Arrays.asList("chrom", "pos", "sample", "GQ", "LofDGivenAA", "LofDGivenAB", "LofDGivenBB", "concordant", "concordantInt", EVAL_NAME, COMP_NAME);
        out.println(Utils.join("\t", fields));
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        if ( tracker.getNBoundRodTracks() > 0 ) {
            for ( String sample : samples ) {
                Genotype evalGT = getGenotype(tracker, ref, sample, EVAL_NAME);
                Genotype compGT = getGenotype(tracker, ref, sample, COMP_NAME);

                if ( evalGT == null && compGT == null && ! keepDoubleMissing )
                    continue;

                if ( (evalGT == null || evalGT.isNoCall()) && ! keepEvalNoCall )
                    continue;

                if ( (compGT == null || compGT.isNoCall()) && ! keepCompNoCall )
                    continue;

                out.printf("%s\t%s\t%s\t", ref.getLocus().getContig(), ref.getLocus().getStart(), sample);
                String evalGQ = evalGT == null ? "NA" : String.format("%d", (int)(10*evalGT.getNegLog10PError()));

                String LofDGivenAA = "-1", LofDGivenAB = "-1", LofDGivenBB = "-1";
                if ( evalGT != null ) {
                    double[] pls = evalGT.getLikelihoods().getAsVector();
                    if ( pls != null ) { // not missing
                        LofDGivenAA = String.format("%.0f", pls[0]);
                        LofDGivenAB = String.format("%.0f", pls[1]);
                        LofDGivenBB = String.format("%.0f", pls[2]);
                    }
                }

                String concordance = evalGT == null || compGT == null ? "NA" : String.format("%s", evalGT.getType() == compGT.getType());
                String concordanceInt = Integer.toString(evalGT == null || compGT == null ? -1 : (evalGT.getType() == compGT.getType() ? 1 : 0));
                String evalType = evalGT == null ? "MISSING" : evalGT.getType().toString();
                String compType = compGT == null ? "MISSING" : compGT.getType().toString();
                out.println(Utils.join("\t", Arrays.asList(evalGQ, LofDGivenAA, LofDGivenAB, LofDGivenBB, concordance, concordanceInt, evalType, compType)));
            }
        }

        return 1;
    }

    private Genotype getGenotype(RefMetaDataTracker tracker, ReferenceContext ref, String sample, String rod) {
        for ( VariantContext vc : tracker.getVariantContexts(ref, rod, null, ref.getLocus(), true, false) ) {
            if ( vc.isNotFiltered() && vc.hasGenotype(sample) )
                return vc.getGenotype(sample);
            else
                return null;
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    @Override
    public void onTraversalDone(Integer sum) {}
}
