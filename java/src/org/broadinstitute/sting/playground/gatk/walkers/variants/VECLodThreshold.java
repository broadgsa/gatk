package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.LocusContext;

public class VECLodThreshold implements VariantExclusionCriterion {
    private double lodThreshold = 5.0;

    private boolean exclude;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            lodThreshold = Double.valueOf(arguments);
        }
    }

    public void compute(char ref, LocusContext context, rodVariants variant) {
        exclude = (variant.getLodBtr() < lodThreshold);
    }

    public boolean useZeroQualityReads() { return false; }

    public boolean exclude(char ref, LocusContext context, rodVariants variant) {
        return (variant.getLodBtr() < lodThreshold);
    }

    public boolean isExcludable() {
        return exclude;
    }

    public String getStudyHeader() {
        return "";
    }

    public String getStudyInfo() {
        return "";
    }

}
