package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;

public class VECNull implements VariantExclusionCriterion {
    public void initialize(String arguments) {
    }

    public void compute(char ref, AlignmentContext context, rodVariants variant) {
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return 1.0;
    }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }

    public boolean useZeroQualityReads() {
        return false;
    }
}
