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

import com.google.java.contract.Ensures;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.Collection;
import java.util.LinkedList;

/**
 * A band pass filtering version of the activity profile
 *
 * Applies a band pass filter with a Gaussian kernel to the input state probabilities to smooth
 * them out of an interval
 *
 * @author Mark DePristo
 * @since 2011
 */
public class BandPassActivityProfile extends ActivityProfile {
    public static final int DEFAULT_FILTER_SIZE = 80;
    public static final double DEFAULT_SIGMA = 55.0;

    private final int filterSize;
    private final double sigma;
    private final double[] GaussianKernel;

    /**
     * Create a band pass activity profile with the default band size
     * @param parser our genome loc parser
     */
    public BandPassActivityProfile(final GenomeLocParser parser) {
        this(parser, DEFAULT_FILTER_SIZE, DEFAULT_SIGMA);
    }

    /**
     * Create an activity profile that implements a band pass filter on the states
     * @param parser our genome loc parser
     * @param filterSize the size (in bp) of the band pass filter.  The filter size is the number of bp to each
     *                   side that are included in the band.  So a filter size of 1 implies that the actual band
     *                   is 3 bp, 1 for the center site and 1 on each size.  2 => 5, etc.
     */
    public BandPassActivityProfile(final GenomeLocParser parser, final int filterSize, final double sigma) {
        super(parser);

        if ( filterSize < 0 ) throw new IllegalArgumentException("Filter size must be greater than or equal to 0 but got " + filterSize);
        if ( sigma < 0 ) throw new IllegalArgumentException("Sigma must be greater than or equal to 0 but got " + sigma);

        // setup the Gaussian kernel for the band pass filter
        this.filterSize = filterSize;
        this.sigma = sigma;
        final double[] kernel = new double[getBandSize()];
        for( int iii = 0; iii < 2* filterSize + 1; iii++ ) {
            kernel[iii] = MathUtils.NormalDistribution(filterSize, sigma, iii);
        }
        this.GaussianKernel = MathUtils.normalizeFromRealSpace(kernel);
    }

    /**
     * Our maximize propagation distance is whatever our parent's is, plus our filter size
     *
     * Stops the profile from interpreting sites that aren't yet fully determined due to
     * propagation of the probabilities.
     *
     * @return the distance in bp we might move our probabilities around for some site i
     */
    @Override
    public int getMaxProbPropagationDistance() {
        return super.getMaxProbPropagationDistance() + filterSize;
    }

    /**
     * Get the size (in bp) of the band pass filter
     * @return a positive integer
     */
    @Ensures("result >= 1")
    public int getBandSize() {
        return 2 * filterSize + 1;
    }

    /**
     * Get the filter size (which is the size of each wing of the band, minus the center point)
     * @return a positive integer
     */
    @Ensures("result >= 0")
    public int getFilteredSize() {
        return filterSize;
    }

    /**
     * Get the Gaussian kernel sigma value
     * @return a positive double
     */
    @Ensures("result >= 0")
    public double getSigma() {
        return sigma;
    }

    /**
     * Get the kernel of this band pass filter.  Do not modify returned result
     * @return the kernel used in this band pass filter
     */
    @Ensures({"result != null", "result.length == getBandSize()"})
    protected double[] getKernel() {
        return GaussianKernel;
    }

    /**
     * Band pass the probabilities in the ActivityProfile, producing a new profile that's band pass filtered
     * @return a new double[] that's the band-pass filtered version of this profile
     */
    @Override
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        final Collection<ActivityProfileState> states = new LinkedList<ActivityProfileState>();

        for ( final ActivityProfileState superState : super.processState(justAddedState) ) {
            for( int jjj = -filterSize; jjj <= filterSize; jjj++ ) {
                final GenomeLoc loc = getLocForOffset(justAddedState.getLoc(), jjj);
                if ( loc != null ) {
                    final double newProb = superState.isActiveProb * GaussianKernel[jjj + filterSize];
                    states.add(new ActivityProfileState(loc, newProb));
                }
            }
        }

        return states;
    }
}
