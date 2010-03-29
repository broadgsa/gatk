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
            // ignore positions at the beginning or end of the overall interval (where there aren't enough records)
            if ( variants[i] == null || variants[i+snpThreshold-1] == null )
                continue;

            // note: not all calls are variant, so we need to be careful.
            // if we don't start with a variant, skip to the next one
            if ( !variants[i].getVariantContext().isVariant() )
                continue;

            // find the nth variant
            GenomeLoc left = variants[i].getVariantContext().getLocation();
            GenomeLoc right = null;
            int snpsSeen = 1;

            int currentIndex = i;
            while ( ++currentIndex < variants.length ) {
                if ( variants[currentIndex].getVariantContext().isVariant() ) {
                    if ( ++snpsSeen == snpThreshold ) {
                        right = variants[currentIndex].getVariantContext().getLocation();
                        break;
                    }
                }
            }

            if ( right != null &&
                 left.getContigIndex() == right.getContigIndex() &&
                 Math.abs(right.getStart() - left.getStart()) <= window )
                return true;
        }
        return false;
    }
}
