/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.VariantContext;

@Analysis(description = "Ti/Tv Variant Evaluator")
public class TiTvVariantEvaluator extends VariantEvaluator implements StandardEval {
    @DataPoint(description = "number of transition loci", format = "%d")
    public long nTi = 0;
    @DataPoint(description = "number of transversion loci", format = "%d")
    public long nTv = 0;
    @DataPoint(description = "the transition to transversion ratio", format = "%.2f")
    public double tiTvRatio = 0.0;
    @DataPoint(description = "number of comp transition sites", format = "%d")
    public long nTiInComp = 0;
    @DataPoint(description = "number of comp transversion sites", format = "%d")
    public long nTvInComp = 0;
    @DataPoint(description = "the transition to transversion ratio for comp sites", format = "%.2f")
    public double TiTvRatioStandard = 0.0;
    @DataPoint(description = "number of derived transition loci", format = "%d")
    public long nTiDerived = 0;
    @DataPoint(description = "number of derived transversion loci", format = "%d")
    public long nTvDerived = 0;
    @DataPoint(description = "the derived transition to transversion ratio", format = "%.2f")
    public double tiTvDerivedRatio = 0.0;

    public int getComparisonOrder() {
        return 2;   // we only need to see each eval track
    }

    public void updateTiTv(VariantContext vc, boolean updateStandard) {
        if (vc != null && vc.isSNP() && vc.isBiallelic() && vc.isPolymorphicInSamples()) {
            if ( GATKVariantContextUtils.isTransition(vc)) {
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

    @Override
    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (eval != null)
            updateTiTv(eval, false);
        if (comp != null)
            updateTiTv(comp, true);
    }

    @Override
    public void finalizeEvaluation() {
        // the ti/tv ratio needs to be set (it's not calculated per-variant).
        this.tiTvRatio = rate(nTi,nTv);
        this.tiTvDerivedRatio = rate(nTiDerived,nTvDerived);
        this.TiTvRatioStandard = rate(nTiInComp, nTvInComp);
    }
}
