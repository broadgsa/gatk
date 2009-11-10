package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;

import java.util.HashMap;

public class VECLodThreshold implements VariantExclusionCriterion {
    private double lodThreshold = 5.0;
    private double lod;
    private boolean exclude;

    public void initialize(HashMap<String,String> args) {
        if ( args.get("lod") != null )
            lodThreshold = Double.valueOf(args.get("lod"));
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        lod = context.getVariant().getLodBtr();
        exclude = lod < lodThreshold;
    }

    public boolean useZeroQualityReads() { return false; }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "LodThreshold("+lodThreshold+")\tlod";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + lod;
    }

    public String getVCFFilterString() {
        return "LOD";
    }

    public String getScoreString() {
        return String.format("%.2f",lod);
    }
}
