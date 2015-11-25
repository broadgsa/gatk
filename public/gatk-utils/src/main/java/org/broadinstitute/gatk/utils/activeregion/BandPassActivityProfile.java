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

package org.broadinstitute.gatk.utils.activeregion;

import com.google.java.contract.Ensures;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.MathUtils;

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
    public static final int MAX_FILTER_SIZE = 50;
    private final static double MIN_PROB_TO_KEEP_IN_FILTER = 1e-5;
    public static final double DEFAULT_SIGMA = 17.0;

    private final int filterSize;
    private final double sigma;
    private final double[] GaussianKernel;

    /**
     * Create a new BandPassActivityProfile with default sigma and filter sizes
     *
     * @see #BandPassActivityProfile(org.broadinstitute.gatk.utils.GenomeLocParser, org.broadinstitute.gatk.utils.GenomeLocSortedSet, int, double, int, double, boolean)
     */
    public BandPassActivityProfile(final GenomeLocParser parser, final GenomeLocSortedSet restrictToIntervals,
                                   final int maxProbPropagationDistance, final double activeProbThreshold) {
        this(parser, restrictToIntervals, maxProbPropagationDistance, activeProbThreshold, MAX_FILTER_SIZE, DEFAULT_SIGMA);
    }

    /**
     * @see #BandPassActivityProfile(org.broadinstitute.gatk.utils.GenomeLocParser, org.broadinstitute.gatk.utils.GenomeLocSortedSet, int, double, int, double, boolean)
     *
     * sets adaptiveFilterSize to true
     */
    public BandPassActivityProfile(final GenomeLocParser parser, final GenomeLocSortedSet restrictToIntervals,
                                   final int maxProbPropagationDistance, final double activeProbThreshold,
                                   final int maxFilterSize, final double sigma) {
        this(parser, restrictToIntervals, maxProbPropagationDistance, activeProbThreshold, maxFilterSize, sigma, true);
    }

    /**
     * Create an activity profile that implements a band pass filter on the states
     *
     * @param parser our genome loc parser
     * @param restrictToIntervals only include states that are within these intervals, if not null
     * @param maxProbPropagationDistance region probability propagation distance beyond it's maximum size
     * @param activeProbThreshold  threshold for the probability of a profile state being active
     * @param maxFilterSize the maximum size of the band pass filter we are allowed to create, regardless of sigma
     * @param sigma the variance of the Gaussian kernel for this band pass filter
     * @param adaptiveFilterSize if true, use the kernel itself to determine the best filter size
     */
    public BandPassActivityProfile(final GenomeLocParser parser, final GenomeLocSortedSet restrictToIntervals, final int maxProbPropagationDistance,
                                   final double activeProbThreshold, final int maxFilterSize, final double sigma, final boolean adaptiveFilterSize) {
        super(parser, maxProbPropagationDistance, activeProbThreshold, restrictToIntervals);

        if ( sigma < 0 ) throw new IllegalArgumentException("Sigma must be greater than or equal to 0 but got " + sigma);

        // setup the Gaussian kernel for the band pass filter
        this.sigma = sigma;
        final double[] fullKernel = makeKernel(maxFilterSize, sigma);
        this.filterSize = adaptiveFilterSize ? determineFilterSize(fullKernel, MIN_PROB_TO_KEEP_IN_FILTER) : maxFilterSize;
        this.GaussianKernel = makeKernel(this.filterSize, sigma);
    }

    protected static int determineFilterSize(final double[] kernel, final double minProbToKeepInFilter) {
        final int middle = (kernel.length - 1) / 2;
        int filterEnd = middle;
        while ( filterEnd > 0 ) {
            if ( kernel[filterEnd - 1] < minProbToKeepInFilter ) {
                break;
            }
            filterEnd--;
        }
        return middle - filterEnd;
    }

    protected static double[] makeKernel(final int filterSize, final double sigma) {
        final int bandSize = 2 * filterSize + 1;
        final double[] kernel = new double[bandSize];
        for( int iii = 0; iii < bandSize; iii++ ) {
            kernel[iii] = MathUtils.normalDistribution(filterSize, sigma, iii);
        }
        return MathUtils.normalizeFromRealSpace(kernel);
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
            if ( superState.isActiveProb > 0.0 ) {
                for( int jjj = -filterSize; jjj <= filterSize; jjj++ ) {
                    final GenomeLoc loc = getLocForOffset(justAddedState.getLoc(), jjj);
                    if ( loc != null ) {
                        final double newProb = superState.isActiveProb * GaussianKernel[jjj + filterSize];
                        states.add(new ActivityProfileState(loc, newProb));
                    }
                }
            } else {
                states.add(justAddedState);
            }
        }

        return states;
    }
}
