/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.activeregion;

import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 *
 *
 * @author Mark DePristo
 * @since 2011
 */
public class BandPassActivityProfile extends ActivityProfile {
    private static final int FILTER_SIZE = 80;
    private static final double[] GaussianKernel;

    static {
        GaussianKernel = new double[2*FILTER_SIZE + 1];
        for( int iii = 0; iii < 2*FILTER_SIZE + 1; iii++ ) {
            GaussianKernel[iii] = MathUtils.NormalDistribution(FILTER_SIZE, 55.0, iii);
        }
    }

    public BandPassActivityProfile(final GenomeLocParser parser) {
        this(parser, new ArrayList<ActivityProfileState>(), null);
    }

    public BandPassActivityProfile(final GenomeLocParser parser, final List<ActivityProfileState> isActiveList, final GenomeLoc regionStartLoc) {
        super(parser, isActiveList, regionStartLoc);
    }

    @Override
    protected ActivityProfile createDerivedProfile(List<ActivityProfileState> isActiveList) {
        return new BandPassActivityProfile(parser,  isActiveList, regionStartLoc);
    }

    /**
     * Band pass the probabilities in the ActivityProfile, producing a new profile that's band pass filtered
     * @return a new double[] that's the band-pass filtered version of this profile
     */
    @Override
    public double[] finalizeProbabilities() {
        final double[] activeProbArray = super.finalizeProbabilities();
        final double[] bandPassProbArray = new double[activeProbArray.length];

        // apply the band pass filter for activeProbArray into filteredProbArray
        for( int iii = 0; iii < activeProbArray.length; iii++ ) {
            final double[] kernel = ArrayUtils.subarray(GaussianKernel, Math.max(FILTER_SIZE-iii, 0), Math.min(GaussianKernel.length,FILTER_SIZE + activeProbArray.length - iii));
            final double[] activeProbSubArray = ArrayUtils.subarray(activeProbArray, Math.max(0,iii - FILTER_SIZE), Math.min(activeProbArray.length,iii + FILTER_SIZE + 1));
            bandPassProbArray[iii] = MathUtils.dotProduct(activeProbSubArray, kernel);
        }

        return bandPassProbArray;
    }
}
