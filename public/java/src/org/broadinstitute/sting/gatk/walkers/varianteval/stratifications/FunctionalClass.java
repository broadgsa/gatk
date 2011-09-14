package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.SnpEff;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by low-, moderate-, and high-impact genomic effect using SnpEff annotations produced by VariantAnnotator
 */
public class FunctionalClass extends VariantStratifier {

    public static final String LOW_IMPACT_STATE_NAME =      "low-impact";
    public static final String MODERATE_IMPACT_STATE_NAME = "moderate-impact";
    public static final String HIGH_IMPACT_STATE_NAME =     "high-impact";

    public static final String EFFECT_IMPACT_ATTRIBUTE_KEY = SnpEff.InfoFieldKey.EFF_IMPACT.toString();

    @Override
    public void initialize() {
        states.add("all");
        states.add(LOW_IMPACT_STATE_NAME);
        states.add(MODERATE_IMPACT_STATE_NAME);
        states.add(HIGH_IMPACT_STATE_NAME);
    }


    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("all");

        if ( eval != null && eval.isVariant() && eval.hasAttribute(EFFECT_IMPACT_ATTRIBUTE_KEY) ) {
            String effectImpact = eval.getAttributeAsString(EFFECT_IMPACT_ATTRIBUTE_KEY);

            if ( effectImpact.equals(SnpEff.EffectImpact.LOW.toString()) ) {
                relevantStates.add(LOW_IMPACT_STATE_NAME);
            }
            else if ( effectImpact.equals(SnpEff.EffectImpact.MODERATE.toString()) ) {
                relevantStates.add(MODERATE_IMPACT_STATE_NAME);
            }
            else if ( effectImpact.equals(SnpEff.EffectImpact.HIGH.toString()) ) {
                relevantStates.add(HIGH_IMPACT_STATE_NAME);
            }
        }

        return relevantStates;
    }
}
