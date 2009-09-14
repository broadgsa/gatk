package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;

import java.util.HashMap;


/**
 * Created by IntelliJ IDEA.
 * User: michaelmelgar
 * Date: Jun 22, 2009
 * Time: 6:04:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class VECDepthOfCoverage implements VariantExclusionCriterion {
    private int maximum = 200;

    private boolean exclude = false;
    private int depth;

    public void initialize(HashMap<String,String> args) {
        if ( args.get("max") != null )
            maximum = Integer.valueOf(args.get("max"));
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        depth = context.getAlignmentContext(useZeroQualityReads()).getReads().size();
        exclude = depth > maximum;
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

//    public boolean isExcludable() {
//        return exclude;
//    }

    public String getStudyHeader() {
        return "DepthOfCoverage("+maximum+")\tdepth";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + depth;
    }

    public String getVCFFilterString() {
        return "DoC";
    }

    public boolean useZeroQualityReads() { return false; }
}
