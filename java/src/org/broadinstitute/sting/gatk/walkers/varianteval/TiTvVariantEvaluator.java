package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;

@Analysis(name = "Ti/Tv Variant Evaluator", description = "Ti/Tv Variant Evaluator")
public class TiTvVariantEvaluator extends VariantEvaluator {

    @DataPoint(name = "ti_count", description = "number of transition loci")
    long nTi = 0;
    @DataPoint(name = "tv_count", description = "number of transversion loci")
    long nTv = 0;
    @DataPoint(name = "ti/tv ratio", description = "the transition to transversion ratio")
    double tiTvRatio = 0.0;    
    @DataPoint(name = "ti_count_std", description = "number of standard transition sites")
    long nTiInStd = 0;
    @DataPoint(name = "tv_count_std", description = "number of standard transversion sites")
    long nTvInStd = 0;
    @DataPoint(name = "ti/tv ratio standard", description = "the transition to transversion ratio")
    double TiTvRatioStandard = 0.0;

    public TiTvVariantEvaluator(VariantEvalWalker parent) {
        super(parent);
    }

    public boolean enabled() {
        return true;
    }

    public String getName() {
        return "titv";
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if (vc != null && vc.isSNP() && vc.isBiallelic()) {
            if (vc.isTransition()) {
                if (updateStandard) nTiInStd++;
                else nTi++;
            } else {                                
                if (updateStandard) nTvInStd++;
                else nTv++;
            }
        }
    }

    public String update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc1 != null) updateTiTv(vc1, false);
        if (vc2 != null) updateTiTv(vc2, true);

        return null; // we don't capture any intersting sites
    }

    @Override
    public void finalizeEvaluation() {
        // the ti/tv ratio needs to be set (it's not calculated per-variant).
        this.tiTvRatio = rate(nTi,nTv);
        this.TiTvRatioStandard = rate(nTiInStd,nTvInStd);
    }
}