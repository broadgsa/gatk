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
import org.broadinstitute.gatk.tools.walkers.varianteval.util.SortableJexlVCMatchExp;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;

/**
 * Stratifies the eval RODs by user-supplied JEXL expressions
 *
 * See http://gatkforums.broadinstitute.org/discussion/1255/what-are-jexl-expressions-and-how-can-i-use-them-with-the-gatk for more details
 */
public class JexlExpression extends VariantStratifier implements StandardStratification {
    // needs to know the jexl expressions
    private Set<SortableJexlVCMatchExp> jexlExpressions;

    @Override
    public void initialize() {
        jexlExpressions = getVariantEvalWalker().getJexlExpressions();

        states.add("none");
        for ( SortableJexlVCMatchExp jexlExpression : jexlExpressions ) {
            states.add(jexlExpression.name);
        }
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        ArrayList<Object> relevantStates = new ArrayList<Object>();
        relevantStates.add("none");

        for ( SortableJexlVCMatchExp jexlExpression : jexlExpressions ) {
            if (eval != null && VariantContextUtils.match(eval, jexlExpression)) {
                relevantStates.add(jexlExpression.name);
            }
        }

        return relevantStates;
    }
}
