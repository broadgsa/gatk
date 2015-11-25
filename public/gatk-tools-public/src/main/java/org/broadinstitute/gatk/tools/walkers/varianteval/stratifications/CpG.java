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

import java.util.ArrayList;
import java.util.List;

/**
 * CpG is a stratification module for VariantEval that divides the input data by within/not within a CpG site
 *
 * <p>
 * It is a three-state stratification:
 * <ul>
 *     <li>The locus is a CpG site ("CpG")
 *     <li>The locus is not a CpG site ("non_CpG")
 *     <li>The locus is either a CpG or not a CpG site ("all")
 * </ul>
 * A CpG site is defined as a site where the reference base at a locus is a C and the adjacent reference base in the 3' direction is a G.
 */
public class CpG extends VariantStratifier {
    @Override
    public void initialize() {
        states.add("all");
        states.add("CpG");
        states.add("non_CpG");
    }

    @Override
    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName, String FamilyName) {
        boolean isCpG = false;
        if (ref != null && ref.getBases() != null) {
            String fwRefBases = new String(ref.getBases());

            //String leftFlank = fwRefBases.substring((fwRefBases.length()/2) - 1, (fwRefBases.length()/2) + 1);
            String rightFlank = fwRefBases.substring((fwRefBases.length()/2), (fwRefBases.length()/2) + 2);

            //if (leftFlank.equalsIgnoreCase("CG") || leftFlank.equalsIgnoreCase("GC") || rightFlank.equalsIgnoreCase("CG") || rightFlank.equalsIgnoreCase("GC")) {
            if (rightFlank.equalsIgnoreCase("CG")) {
                isCpG = true;
            }
        }

        ArrayList<Object> relevantStates = new ArrayList<Object>(2);
        relevantStates.add("all");
        relevantStates.add(isCpG ? "CpG" : "non_CpG");

        return relevantStates;
    }
}
