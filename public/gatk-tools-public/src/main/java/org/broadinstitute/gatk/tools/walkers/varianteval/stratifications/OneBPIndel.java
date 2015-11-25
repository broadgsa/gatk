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

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the eval RODs into sites where the indel is 1 bp in length and those where the event is 2+.
 * all non indel events go into all bins, so that SNP counts can be used as contrasts in eval modules.
 */
public class OneBPIndel extends VariantStratifier {
    private final static List<Object> ALL = Arrays.asList((Object)"all", (Object)"one.bp", (Object)"two.plus.bp");
    private final static List<Object> ONE_BP = Arrays.asList((Object)"all", (Object)"one.bp");
    private final static List<Object> TWO_PLUS_BP = Arrays.asList((Object)"all", (Object)"two.plus.bp");

    @Override
    public void initialize() {
        states.addAll(ALL);
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (eval != null && eval.isIndel()) {
            for ( int l : eval.getIndelLengths() )
                if ( Math.abs(l) > 1 )
                    return TWO_PLUS_BP; // someone is too long
            return ONE_BP; // all lengths are one
        } else
            return ALL;
    }
}
