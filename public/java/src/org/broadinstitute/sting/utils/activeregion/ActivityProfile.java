/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.activeregion;

import com.google.java.contract.Requires;
import org.apache.commons.lang.ArrayUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Class holding information about per-base activity scores for the
 * active region traversal
 *
 * @author Mark DePristo
 * @since Date created
 */
public class ActivityProfile {
    final GenomeLocParser parser;
    final boolean presetRegions;
    GenomeLoc regionStartLoc = null;
    GenomeLoc regionStopLoc = null;
    final List<ActivityProfileResult> isActiveList;
    private static final int FILTER_SIZE = 80;
    private static final double[] GaussianKernel;

    static {
        GaussianKernel = new double[2*FILTER_SIZE + 1];
        for( int iii = 0; iii < 2*FILTER_SIZE + 1; iii++ ) {
            GaussianKernel[iii] = MathUtils.NormalDistribution(FILTER_SIZE, 55.0, iii);
        }
    }

    // todo -- add upfront the start and stop of the intervals
    // todo -- check that no regions are unexpectedly missing
    // todo -- add unit tests
    // TODO -- own preset regions
    public ActivityProfile(final GenomeLocParser parser, final boolean presetRegions) {
        this(parser, presetRegions, new ArrayList<ActivityProfileResult>(), null);
    }

    protected ActivityProfile(final GenomeLocParser parser, final boolean presetRegions, final List<ActivityProfileResult> isActiveList, final GenomeLoc regionStartLoc) {
        this.parser = parser;
        this.presetRegions = presetRegions;
        this.isActiveList = isActiveList;
        this.regionStartLoc = regionStartLoc;
    }

    @Override
    public String toString() {
        return "ActivityProfile{" +
                "start=" + regionStartLoc +
                ", stop=" + regionStopLoc +
                '}';
    }

    /**
     * Add the next ActivityProfileResult to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param result a well-formed ActivityProfileResult result to incorporate into this profile
     */
    @Requires("result != null")
    public void add(final ActivityProfileResult result) {
        final GenomeLoc loc = result.getLoc();

        if ( regionStartLoc == null ) {
            regionStartLoc = loc;
            regionStopLoc = loc;
        } else {
            if ( regionStopLoc.getStart() != loc.getStart() - 1 )
                throw new IllegalArgumentException("Bad add call to ActivityProfile: loc " + loc + " not immediate after last loc " + regionStopLoc );
            regionStopLoc = loc;
        }

        isActiveList.add(result);
    }

    public int size() {
        return isActiveList.size();
    }

    public boolean isEmpty() {
        return isActiveList.isEmpty();
    }

    public boolean hasPresetRegions() {
        return presetRegions;
    }

    /**
     * Band pass this ActivityProfile, producing a new profile that's band pass filtered
     * @return a new ActivityProfile that's the band-pass filtered version of this profile
     */
    public ActivityProfile bandPassFilter() {
        final double[] activeProbArray = new double[isActiveList.size()];
        int iii = 0;
        for( final ActivityProfileResult result : isActiveList ) {
            activeProbArray[iii++] = result.isActiveProb;
        }
        iii = 0;
        for( final ActivityProfileResult result : isActiveList ) {
            if( result.resultState.equals(ActivityProfileResult.ActivityProfileResultState.HIGH_QUALITY_SOFT_CLIPS) ) { // special code to deal with the problem that high quality soft clipped bases aren't added to pileups
                final int numHQClips = result.resultValue.intValue();
                for( int jjj = Math.max(0, iii - numHQClips); jjj < Math.min(activeProbArray.length, iii+numHQClips); jjj++ ) {
                    activeProbArray[jjj] = Math.max(activeProbArray[jjj], activeProbArray[iii]);
                }
            }
            iii++;
        }

        final double[] filteredProbArray;
        if( !presetRegions ) {
            // if we aren't using preset regions, actually apply the band pass filter for activeProbArray into filteredProbArray
            filteredProbArray = new double[activeProbArray.length];
            for( iii = 0; iii < activeProbArray.length; iii++ ) {
                final double[] kernel = ArrayUtils.subarray(GaussianKernel, Math.max(FILTER_SIZE-iii, 0), Math.min(GaussianKernel.length,FILTER_SIZE + activeProbArray.length - iii));
                final double[] activeProbSubArray = ArrayUtils.subarray(activeProbArray, Math.max(0,iii - FILTER_SIZE), Math.min(activeProbArray.length,iii + FILTER_SIZE + 1));
                filteredProbArray[iii] = MathUtils.dotProduct(activeProbSubArray, kernel);
            }
        } else {
            // otherwise we simply use the activeProbArray directly
            filteredProbArray = activeProbArray;
        }

        iii = 0;
        for( final double prob : filteredProbArray ) {
            final ActivityProfileResult result = isActiveList.get(iii++);
            result.isActiveProb = prob;
            result.resultState = ActivityProfileResult.ActivityProfileResultState.NONE;
            result.resultValue = null;
        }

        return new ActivityProfile(parser, presetRegions, isActiveList, regionStartLoc);
    }

    /**
     * Partition this profile into active regions
     * @param activeRegionExtension the amount of margin overlap in the active region
     * @return the list of active regions
     */
    public List<ActiveRegion> createActiveRegions( final int activeRegionExtension, final int maxRegionSize ) {
        final double ACTIVE_PROB_THRESHOLD = 0.002; // TODO: needs to be set-able by the walker author
        final ArrayList<ActiveRegion> returnList = new ArrayList<ActiveRegion>();

        if( isActiveList.size() == 0 ) {
            // no elements in the active list, just return an empty one
            return Collections.emptyList();
        } else if( isActiveList.size() == 1 ) {
            // there's a single element, it's either active or inactive
            boolean isActive = isActiveList.get(0).isActiveProb > ACTIVE_PROB_THRESHOLD;
            returnList.addAll(createActiveRegion(isActive, 0, 0, activeRegionExtension, maxRegionSize));
        } else {
            // there are 2+ elements, divide these up into regions
            boolean isActive = isActiveList.get(0).isActiveProb > ACTIVE_PROB_THRESHOLD;
            int curStart = 0;
            for(int iii = 1; iii < isActiveList.size(); iii++ ) {
                final boolean thisStatus = isActiveList.get(iii).isActiveProb > ACTIVE_PROB_THRESHOLD;
                if( isActive != thisStatus ) {
                    returnList.addAll(createActiveRegion(isActive, curStart, iii - 1, activeRegionExtension, maxRegionSize));
                    isActive = thisStatus;
                    curStart = iii;
                }
            }
            returnList.addAll(createActiveRegion(isActive, curStart, isActiveList.size() - 1, activeRegionExtension, maxRegionSize)); // close out the current active region
        }
        return returnList;
    }

    /**
     * Helper routine to create an active region based on our current start and end offsets
     * @param isActive should the region be active?
     * @param curStart offset (0-based) from the start of this region
     * @param curEnd offset (0-based) from the start of this region
     * @param activeRegionExtension the amount of margin overlap in the active region
     * @return a fully initialized ActiveRegion with the above properties
     */
    private final List<ActiveRegion> createActiveRegion(final boolean isActive, final int curStart, final int curEnd, final int activeRegionExtension, final int maxRegionSize) {
        return createActiveRegion(isActive, curStart, curEnd, activeRegionExtension, maxRegionSize, new ArrayList<ActiveRegion>());
    }

    private final List<ActiveRegion> createActiveRegion(final boolean isActive, final int curStart, final int curEnd, final int activeRegionExtension, final int maxRegionSize, final List<ActiveRegion> returnList) {
        if( !isActive || curEnd - curStart < maxRegionSize ) {
            final GenomeLoc loc = parser.createGenomeLoc(regionStartLoc.getContig(), regionStartLoc.getStart() + curStart, regionStartLoc.getStart() + curEnd);
            returnList.add(new ActiveRegion(loc, isActive, parser, activeRegionExtension));
            return returnList;
        }
        // find the best place to break up the large active region
        Double minProb = Double.MAX_VALUE;
        int cutPoint = -1;

        final int size = curEnd - curStart + 1;
        for( int iii = curStart + (int)(size*0.15); iii < curEnd - (int)(size*0.15); iii++ ) {
            if( isActiveList.get(iii).isActiveProb < minProb ) { minProb = isActiveList.get(iii).isActiveProb; cutPoint = iii; }
        }
        final List<ActiveRegion> leftList = createActiveRegion(isActive, curStart, cutPoint, activeRegionExtension, maxRegionSize, new ArrayList<ActiveRegion>());
        final List<ActiveRegion> rightList = createActiveRegion(isActive, cutPoint+1, curEnd, activeRegionExtension, maxRegionSize, new ArrayList<ActiveRegion>());
        returnList.addAll( leftList );
        returnList.addAll( rightList );
        return returnList;
    }
}
