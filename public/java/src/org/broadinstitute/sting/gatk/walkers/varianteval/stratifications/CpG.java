package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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

    public List<String> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
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

        ArrayList<String> relevantStates = new ArrayList<String>();
        relevantStates.add("all");
        relevantStates.add(isCpG ? "CpG" : "non_CpG");

        return relevantStates;
    }
}
