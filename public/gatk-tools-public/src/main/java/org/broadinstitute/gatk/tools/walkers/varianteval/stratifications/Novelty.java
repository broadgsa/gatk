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

import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.*;

/**
 * Stratifies by whether a site in in the list of known RODs (e.g., dbsnp by default)
 */
public class Novelty extends VariantStratifier implements StandardStratification {
    // needs the variant contexts and known names
    private List<RodBinding<VariantContext>> knowns;

    private final static List<Object> KNOWN_STATES = Arrays.asList((Object)"all", (Object)"known");
    private final static List<Object> NOVEL_STATES = Arrays.asList((Object)"all", (Object)"novel");

    @Override
    public void initialize() {
        states.addAll(Arrays.asList("all", "known", "novel"));
        knowns = getVariantEvalWalker().getKnowns();
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        if (tracker != null && eval != null) {
            final Collection<VariantContext> knownComps = tracker.getValues(knowns, ref.getLocus());
            for ( final VariantContext c : knownComps ) {
                // loop over sites, looking for something that matches the type eval
                if ( eval.getType() == c.getType() || eval.getType() == VariantContext.Type.NO_VARIATION ) {
                    return KNOWN_STATES;
                }
            }
        } 
        
        return NOVEL_STATES;
    }
}