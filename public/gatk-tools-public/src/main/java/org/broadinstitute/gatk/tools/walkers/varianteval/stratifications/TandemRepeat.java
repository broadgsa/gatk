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
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.List;

/**
 * Stratifies the eval RODs into sites that are tandem repeats
 */
public class TandemRepeat extends VariantStratifier {
    private final static List<Object> JUST_ALL = Arrays.asList((Object)"all");
    private final static List<Object> ALL = Arrays.asList((Object)"all", (Object)"is.repeat", (Object)"not.repeat");
    private final static List<Object> REPEAT = Arrays.asList((Object)"all", (Object)"is.repeat");
    private final static List<Object> NOT_REPEAT = Arrays.asList((Object)"all", (Object)"not.repeat");

    @Override
    public void initialize() {
        states.addAll(ALL);
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if ( eval == null || ! eval.isIndel() )
            return ALL;
        else if ( GATKVariantContextUtils.isTandemRepeat(eval, ref.getForwardBases()) ) {
            print("REPEAT", eval, ref);
            return REPEAT;
        } else {
            print("NOT A REPEAT", eval, ref);
            return NOT_REPEAT;
        }
    }
    
    private final void print(String prefix, VariantContext eval, ReferenceContext ref) {
//        String alleles = ParsingUtils.sortList(eval.getAlleles()).toString();
//        this.getVariantEvalWalker().getLogger().info(prefix + ": " + "pos=" + eval.getStart() + " alleles=" + alleles + " ref=" + new String(ref.getForwardBases()));
    }
}
