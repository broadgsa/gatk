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
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval RODs by the indel size
 *
 * Indel sizes are stratified from sizes -100 to +100. Sizes greater than this are lumped in the +/- 100 bin
 * This stratification ignores multi-allelic indels (whose size is not defined uniquely)
 */
public class IndelSize extends VariantStratifier {
    static final int MAX_INDEL_SIZE = 100;

    @Override
    public void initialize() {
        for( int a=-MAX_INDEL_SIZE; a <=MAX_INDEL_SIZE; a++ ) {
            states.add(a);
        }
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null && eval.isIndel() && eval.isBiallelic()) {
            try {
                int eventLength = 0;
                if ( eval.isSimpleInsertion() ) {
                    eventLength = eval.getAlternateAllele(0).length();
                } else if ( eval.isSimpleDeletion() ) {
                    eventLength = -eval.getReference().length();
                }

                if (eventLength > MAX_INDEL_SIZE)
                    eventLength = MAX_INDEL_SIZE;
                else if (eventLength < -MAX_INDEL_SIZE)
                    eventLength = -MAX_INDEL_SIZE;

                return Collections.singletonList((Object)eventLength);
            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }
    @Override
    public String getFormat() {
        return "%d";
    }
}
