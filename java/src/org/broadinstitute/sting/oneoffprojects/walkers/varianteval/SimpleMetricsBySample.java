package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.SampleDataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvaluatorBySample;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;

import java.util.ArrayList;
import java.util.List;

/**
 * Extends the Per-sample variant evaluator class and returns, for each sample, the number of variants, the Ti/Tv, and
 * the comp overlap. It does this only on sites where the sample is identified as hom var, or het.
 */
@Analysis(name = "Simple Metrics by Sample", description = "Variant counts, Ti/Tv, comp overlap; per sample")
public class SimpleMetricsBySample extends VariantEvaluatorBySample {
    public SimpleMetricsBySample(VariantEvalWalker parent) { super(parent); }

    public List<SampleDataPoint> getDataPoints() {
        List<SampleDataPoint> points = new ArrayList(3);
        points.add(new CountSNPsSample());
        points.add(new TiTvRatioSample());
        points.add(new CompOverlapSample());

        return points;
    }

    public String getTableName() {
        return "SimpleMetricsBySample";
    }

    public String getName() {
        return "SimpleMetricsBySample";
    }

    public int getComparisonOrder() { return 2; }

    public boolean includeGenotype(Genotype g) {
        return (g.isHet() || g.isHomVar()) && ! g.isFiltered();
    }

    public boolean enabled() {
        return true;
    }

}

class CountSNPsSample extends SampleDataPoint {
    int numVariants = 0;

    public CountSNPsSample() {
        super("CountVariants");
    }

    public void update2(VariantContext vc, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc != null && vc.isSNP() ) {
            numVariants++;
        }
    }

    public String toString() {
        return String.format("%d",numVariants);
    }
}

class TiTvRatioSample extends SampleDataPoint {
    int nTi = 0;
    int nTv = 0;

    public TiTvRatioSample() {
        super("TiTvRatio");
    }

    public void update2(VariantContext vc, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc != null && vc.isSNP() && vc.isBiallelic() ) {
            if ( VariantContextUtils.isTransition(vc) ) {
                nTi++;
            } else {
                nTv++;
            }
        }
    }

    public String toString() {
        return String.format("%.2f", ( ((double) nTi )/ nTv));
    }
}

class CompOverlapSample extends SampleDataPoint {
    int nOverlap = 0;

    public CompOverlapSample() {
        super("CompOverlap");
    }

    public void update2(VariantContext eval, VariantContext comp,RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        boolean compIsGood = comp != null && comp.isNotFiltered() && comp.isSNP() ;
        boolean evalIsGood = eval != null && eval.isSNP();
        if ( compIsGood && evalIsGood ) {
            nOverlap++;
        }
    }

    public String toString() {
        return String.format("%d",nOverlap);
    }
}
