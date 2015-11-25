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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications;

import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.SnpEff;
import htsjdk.variant.variantcontext.VariantContext;

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


    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        ArrayList<Object> relevantStates = new ArrayList<Object>();

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
