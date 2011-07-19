package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

public class ClusteredSnps {
    private GenomeLocParser genomeLocParser;
    private int window = 10;
    private int snpThreshold = 3;

    public ClusteredSnps(GenomeLocParser genomeLocParser,int snpThreshold, int window) {
        this.genomeLocParser = genomeLocParser;
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

            // note: the documentation tells users we'll blow up if ref calls are present.
            //   if we ever get a windowed rod context that isn't a hack, we can actually allow this...
            if ( !variants[i].getVariantContext().isVariant() )
                throw new UserException.BadInput("The clustered SNPs filter does not work in the presence of non-variant records; see the documentation for more details");

            // find the nth variant
            GenomeLoc left = VariantContextUtils.getLocation(genomeLocParser,variants[i].getVariantContext());
            GenomeLoc right = null;
            int snpsSeen = 1;

            int currentIndex = i;
            while ( ++currentIndex < variants.length ) {
                if ( variants[currentIndex] != null && variants[currentIndex].getVariantContext() != null && variants[currentIndex].getVariantContext().isVariant() ) {
                    if ( ++snpsSeen == snpThreshold ) {
                        right = VariantContextUtils.getLocation(genomeLocParser,variants[currentIndex].getVariantContext());
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
