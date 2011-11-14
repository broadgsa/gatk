/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;

@Analysis(description = "Build 1000 Genome Phase I paper summary of variants table")
public class G1KPhaseITable extends VariantEvaluator {
    // basic counts on various rates found
    @DataPoint(description = "Number of samples")
    public long nSamples = 0;

    @DataPoint(description = "Number of processed loci")
    public long nProcessedLoci = 0;

    @DataPoint(description = "Number of SNPs")
    public long nSNPs = 0;
    @DataPoint(description = "SNP Novelty Rate")
    public String SNPNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SNPs per individual")
    public long nSNPsPerSample = 0;

    @DataPoint(description = "Number of Indels")
    public long nIndels = 0;
    @DataPoint(description = "Indel Novelty Rate")
    public String IndelNoveltyRate = "NA";
    @DataPoint(description = "Mean number of Indels per individual")
    public long nIndelsPerSample = 0;

    @DataPoint(description = "Number of SVs")
    public long nSVs = 0;
    @DataPoint(description = "SV Novelty Rate")
    public String SVNoveltyRate = "NA";
    @DataPoint(description = "Mean number of SVs per individual")
    public long nSVsPerSample = 0;

    Map<VariantContext.Type, Integer> allVariantCounts, knownVariantCounts;
    Map<String, Map<VariantContext.Type, Integer>> countsPerSample;

    private final Map<VariantContext.Type, Integer> makeCounts() {
        Map<VariantContext.Type, Integer> counts = new EnumMap<VariantContext.Type, Integer>(VariantContext.Type.class);
        counts.put(VariantContext.Type.SNP, 0);
        counts.put(VariantContext.Type.INDEL, 0);
        counts.put(VariantContext.Type.SYMBOLIC, 0);
        return counts;
    }

    public void initialize(VariantEvalWalker walker) {
        countsPerSample = new HashMap<String, Map<VariantContext.Type, Integer>>();
        nSamples = walker.getSampleNamesForEvaluation().size();

        for ( String sample : walker.getSampleNamesForEvaluation() ) {
            countsPerSample.put(sample, makeCounts());
        }

        allVariantCounts = makeCounts();
        knownVariantCounts = makeCounts();
    }

    @Override public boolean enabled() { return true; }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        nProcessedLoci += context.getSkippedBases() + (ref == null ? 0 : 1);
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || eval.isMonomorphic() ) return null;

        switch (eval.getType()) {
            case SNP:
            case INDEL:
            case SYMBOLIC:
                allVariantCounts.put(eval.getType(), allVariantCounts.get(eval.getType()) + 1);
                if ( comp != null )
                    knownVariantCounts.put(eval.getType(), knownVariantCounts.get(eval.getType()) + 1);
                break;
            default:
                throw new UserException.BadInput("Unexpected variant context type: " + eval);
        }

        // count variants per sample
        for (final Genotype g : eval.getGenotypes()) {
            if ( ! g.isNoCall() && ! g.isHomRef() ) {
                int count = countsPerSample.get(g.getSampleName()).get(eval.getType());
                countsPerSample.get(g.getSampleName()).put(eval.getType(), count + 1);
            }
        }

        return null; // we don't capture any interesting sites
    }

    private final int perSampleMean(VariantContext.Type type) {
        long sum = 0;
        for ( Map<VariantContext.Type, Integer> count : countsPerSample.values() ) {
            sum += count.get(type);
        }
        return (int)(Math.round(sum / (1.0 * countsPerSample.size())));
    }

    private final String noveltyRate(VariantContext.Type type) {
        int all = allVariantCounts.get(type);
        int known = knownVariantCounts.get(type);
        int novel = all - known;
        double rate = (novel / (1.0 * all));
        return all == 0 ? "NA" : String.format("%.2f", rate);
    }

    public void finalizeEvaluation() {
        nSNPs = allVariantCounts.get(VariantContext.Type.SNP);
        nIndels = allVariantCounts.get(VariantContext.Type.INDEL);
        nSVs = allVariantCounts.get(VariantContext.Type.SYMBOLIC);

        nSNPsPerSample = perSampleMean(VariantContext.Type.SNP);
        nIndelsPerSample = perSampleMean(VariantContext.Type.INDEL);
        nSVsPerSample = perSampleMean(VariantContext.Type.SYMBOLIC);

        SNPNoveltyRate = noveltyRate(VariantContext.Type.SNP);
        IndelNoveltyRate = noveltyRate(VariantContext.Type.INDEL);
        SVNoveltyRate = noveltyRate(VariantContext.Type.SYMBOLIC);
    }
}