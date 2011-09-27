package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.SnpEff;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Stratifies by nonsense, missense, silent, and all annotations in the input ROD, from the INFO field annotation.
 */
public class FunctionalClass extends VariantStratifier {

    public enum FunctionalType {
        silent,
        missense,
        nonsense
    }


    @Override
    public void initialize() {
        states.add("all");
        for ( FunctionalType type : FunctionalType.values() )
            states.add(type.name());
    }


public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        ArrayList<String> relevantStates = new ArrayList<String>();

        relevantStates.add("all");

        if (eval != null && eval.isVariant()) {
            FunctionalType type = null;

            if (eval.hasAttribute("refseq.functionalClass")) {
                try {
                    type = FunctionalType.valueOf(eval.getAttributeAsString("refseq.functionalClass", null));
                } catch ( Exception e ) {} // don't error out if the type isn't supported
            } else if (eval.hasAttribute("refseq.functionalClass_1")) {
                int annotationId = 1;
                String key;

                do {
                    key = String.format("refseq.functionalClass_%d", annotationId);

                    String newtypeStr = eval.getAttributeAsString(key, null);
                    if ( newtypeStr != null && !newtypeStr.equalsIgnoreCase("null") ) {
                        try {
                            FunctionalType newType = FunctionalType.valueOf(newtypeStr);
                            if ( type == null ||
                                ( type == FunctionalType.silent && newType != FunctionalType.silent ) ||
                                ( type == FunctionalType.missense && newType == FunctionalType.nonsense ) ) {
                                type = newType;
                            }
                        } catch ( Exception e ) {} // don't error out if the type isn't supported
                    }

                    annotationId++;
                } while (eval.hasAttribute(key));

            } else if ( eval.hasAttribute(SnpEff.InfoFieldKey.FUNCTIONAL_CLASS_KEY.getKeyName()) ) {
                try {
                    SnpEff.EffectFunctionalClass snpEffFunctionalClass = SnpEff.EffectFunctionalClass.valueOf(eval.getAttribute(SnpEff.InfoFieldKey.FUNCTIONAL_CLASS_KEY.getKeyName()).toString());
                    if ( snpEffFunctionalClass == SnpEff.EffectFunctionalClass.NONSENSE )
                        type = FunctionalType.nonsense;
                    else if ( snpEffFunctionalClass == SnpEff.EffectFunctionalClass.MISSENSE )
                        type = FunctionalType.missense;
                    else if ( snpEffFunctionalClass == SnpEff.EffectFunctionalClass.SILENT )
                        type = FunctionalType.silent;
                    }
                catch ( Exception e ) {} // don't error out if the type isn't supported
            }

            if ( type != null ) {
                relevantStates.add(type.name());
            }
        }

        return relevantStates;
    }
}
