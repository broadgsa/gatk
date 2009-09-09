package org.broadinstitute.sting.gatk.walkers.filters;

import java.util.HashMap;

public class VECLodThresholdByCoverage implements VariantExclusionCriterion {
    private double slope = 0.46;
    private double depth;
    private double lod;
    private boolean exclude;

    public void initialize(HashMap<String, String> arguments) {
        if (arguments.get("slope") != null) {
            slope = Double.valueOf(arguments.get("slope"));
        }
    }

    public void compute(VariantContextWindow contextWindow) {
        depth = (double) contextWindow.getContext().getVariant().getPileupDepth();
        lod = contextWindow.getContext().getVariant().getLodBtr();

        exclude = (lod < slope*depth);
    }

    public double inclusionProbability() {
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return String.format("LodThresholdByCoverage(%f)\tdepth\tlod", slope);
    }

    public String getStudyInfo() {
        return String.format("%s\t%d\t%f", exclude ? "fail" : "pass", (int) depth, lod);
    }

    public boolean useZeroQualityReads() {
        return false;
    }
}
