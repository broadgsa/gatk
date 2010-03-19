package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.utils.GenomeLoc;

public class ClusteredSnps {
    private int window = 10;
    private int snpThreshold = 3;

    public ClusteredSnps(int snpThreshold, int window) {
        this.window = window;
        this.snpThreshold = snpThreshold;
        if ( window < 1 || snpThreshold < 1 )
            throw new IllegalArgumentException("Window and threshold values need to be positive values");
    }

    public boolean filter(FiltrationContextWindow contextWindow) {

        FiltrationContext[] variants = contextWindow.getWindow(snpThreshold-1, snpThreshold-1);
        for (int i = 0; i < snpThreshold; i++) {
            if ( variants[i] == null || variants[i+snpThreshold-1] == null )
                continue;

            GenomeLoc left = variants[i].getVariantContext().getLocation();
            GenomeLoc right = variants[i+snpThreshold-1].getVariantContext().getLocation();
            if ( left.getContigIndex() == right.getContigIndex() &&
                 Math.abs(right.getStart() - left.getStart()) <= window )
                return true;
        }
        return false;
    }
}
