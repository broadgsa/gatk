package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.LocusContext;

public class VECLodThreshold implements VariantExclusionCriterion {
    private double lodThreshold = 5.0;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            lodThreshold = Double.valueOf(arguments);
        }
    }

    public boolean useZeroQualityReads() { return false; }

    public boolean exclude(char ref, LocusContext context, rodVariants variant) {
        return (variant.getLodBtr() < lodThreshold);
    }
}
