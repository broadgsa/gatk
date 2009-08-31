package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;



public class VECClusteredSnps implements VariantExclusionCriterion {
    private int window = 10;
    private int snpThreshold = 3;
    private boolean exclude;
    private long distance;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            window = Integer.valueOf(argPieces[0]);
            snpThreshold = Integer.valueOf(argPieces[1]);
        }

        if ( window < 2 || snpThreshold < 2 )
            throw new StingException("Window and threshold values need to be >= 2");
    }

    public void compute(VariantContextWindow contextWindow) {
        exclude = false;

        VariantContext[] variants = contextWindow.getWindow(snpThreshold-1, snpThreshold-1);
        for (int i = 0; i < snpThreshold; i++) {
            if ( variants[i] == null || variants[i+snpThreshold-1] == null )
                continue;

            GenomeLoc left = variants[i].getAlignmentContext().getLocation();
            GenomeLoc right = variants[i+snpThreshold-1].getAlignmentContext().getLocation();
            if ( left.getContigIndex() == right.getContigIndex() ) {
                distance = Math.abs(right.getStart() - left.getStart());
                if ( distance <= window ) {
                    exclude = true;
                    return;
                }
            }
        }
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "ClusteredSnps("+window+","+snpThreshold+")\tWindow_size";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + (exclude ? distance : "N/A");
    }

    public boolean useZeroQualityReads() { return false; }
}