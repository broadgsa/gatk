/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.filters;

import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;

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
            GenomeLoc left = GATKVariantContextUtils.getLocation(genomeLocParser, variants[i].getVariantContext());
            GenomeLoc right = null;
            int snpsSeen = 1;

            int currentIndex = i;
            while ( ++currentIndex < variants.length ) {
                if ( variants[currentIndex] != null && variants[currentIndex].getVariantContext() != null && variants[currentIndex].getVariantContext().isVariant() ) {
                    if ( ++snpsSeen == snpThreshold ) {
                        right = GATKVariantContextUtils.getLocation(genomeLocParser, variants[currentIndex].getVariantContext());
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
