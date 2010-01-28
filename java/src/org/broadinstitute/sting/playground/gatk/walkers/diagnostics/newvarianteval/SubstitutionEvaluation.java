package org.broadinstitute.sting.playground.gatk.walkers.diagnostics.newvarianteval;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.playground.utils.NamedTable;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.List;
import java.util.ArrayList;

public class SubstitutionEvaluation implements VariantEvaluation {
    private NamedTable substitutions;

    public SubstitutionEvaluation() {
        substitutions = new NamedTable();
    }

    public void update(RodVCF variant, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (variant != null) {
            List<String> altAlleles = variant.getAlternateAlleleList();
            String altAllele = altAlleles.get(0);
            String refAllele = String.format("%c", BaseUtils.baseIndexToSimpleBase(ref.getBaseIndex()));

            substitutions.increment(refAllele, altAllele);
        }
    }

    public String getName() {
        return "Substitution";
    }

    public List<String> getResult() {
        ArrayList<String> results = new ArrayList<String>();

        int transitions = 0;
        int transversions = 0;
        int unknown = 0;

        for (char rbase : BaseUtils.BASES) {
            String refAllele = String.format("%c", rbase);

            for (char abase : BaseUtils.BASES) {
                String altAllele = String.format("%c", abase);

                if (BaseUtils.isTransition((byte)rbase, (byte)abase)) {
                    transitions += substitutions.get(refAllele, altAllele);
                } else if (BaseUtils.isTransversion((byte)rbase, (byte)abase)) {
                    transversions += substitutions.get(refAllele, altAllele);
                } else {
                    unknown += substitutions.get(refAllele, altAllele);
                }
            }
        }

        results.add(String.format("%n%s", substitutions.toString()));
        results.add(String.format("transitions=%d", transitions));
        results.add(String.format("transversions=%d", transversions));
        results.add(String.format("unknown=%d", unknown));
        results.add(String.format("titv=%f", ((double) transitions)/((double) transversions)));

        return results;
    }
}
