package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

@Analysis(description = "Ti/Tv Variant Evaluator")
public class TiTvVariantEvaluator extends VariantEvaluator implements StandardEval {

    @DataPoint(description = "number of transition loci")
    long nTi = 0;
    @DataPoint(description = "number of transversion loci")
    long nTv = 0;
    @DataPoint(description = "the transition to transversion ratio")
    double tiTvRatio = 0.0;    
    @DataPoint(description = "number of comp transition sites")
    long nTiInComp = 0;
    @DataPoint(description = "number of comp transversion sites")
    long nTvInComp = 0;
    @DataPoint(description = "the transition to transversion ratio for comp sites")
    double TiTvRatioStandard = 0.0;
    @DataPoint(description = "number of derived transition loci")
    long nTiDerived = 0;
    @DataPoint(description = "number of derived transversion loci")
    long nTvDerived = 0;
    @DataPoint(description = "the derived transition to transversion ratio")
    double tiTvDerivedRatio = 0.0;

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if (vc != null && vc.isSNP() && vc.isBiallelic() && vc.isPolymorphic()) {
            if (VariantContextUtils.isTransition(vc)) {
                if (updateStandard) nTiInComp++;
                else nTi++;
            } else {                                
                if (updateStandard) nTvInComp++;
                else nTv++;
            }

            if (vc.hasAttribute("ANCESTRALALLELE")) {
                final String aaStr = vc.getAttributeAsString("ANCESTRALALLELE", "null").toUpperCase();
                if ( ! aaStr.equals(".") ) {
                    switch ( BaseUtils.SNPSubstitutionType(aaStr.getBytes()[0], vc.getAlternateAllele(0).getBases()[0] ) ) {
                        case TRANSITION: nTiDerived++; break;
                        case TRANSVERSION: nTvDerived++; break;
                        default: break;
                    }
                }
            }
        }
    }

    public String update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc1 != null) updateTiTv(vc1, false);
        if (vc2 != null) updateTiTv(vc2, true);

        return null; // we don't capture any interesting sites
    }

    @Override
    public void finalizeEvaluation() {
        // the ti/tv ratio needs to be set (it's not calculated per-variant).
        this.tiTvRatio = rate(nTi,nTv);
        this.tiTvDerivedRatio = rate(nTiDerived,nTvDerived);
        this.TiTvRatioStandard = rate(nTiInComp, nTvInComp);
    }
}
