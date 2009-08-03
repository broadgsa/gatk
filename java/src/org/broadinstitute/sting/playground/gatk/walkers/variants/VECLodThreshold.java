package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.LocusContext;

public class VECLodThreshold implements VariantExclusionCriterion {
    private double lodThreshold = 5.0;
    private double lod;
    private boolean exclude;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            lodThreshold = Double.valueOf(arguments);
        }
    }

    public void compute(char ref, LocusContext context, rodVariants variant) {
        lod = variant.getLodBtr();
        exclude = lod < lodThreshold;
    }

    public boolean useZeroQualityReads() { return false; }

    public boolean isExcludable() {
        return exclude;
    }

    public String getStudyHeader() {
        return "LodThreshold("+lod+")\tlod";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + lod;
    }

}
