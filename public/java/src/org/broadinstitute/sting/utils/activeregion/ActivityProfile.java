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
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

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
    private final static double ACTIVE_PROB_THRESHOLD = 0.002; // TODO: needs to be set-able by the walker author

    protected final List<ActivityProfileState> isActiveList;
    protected final GenomeLocParser parser;

    protected GenomeLoc regionStartLoc = null;
    protected GenomeLoc regionStopLoc = null;

    public ActivityProfile(final GenomeLocParser parser) {
        this(parser, new ArrayList<ActivityProfileState>(), null);
    }

    protected ActivityProfile(final GenomeLocParser parser, final List<ActivityProfileState> isActiveList, final GenomeLoc regionStartLoc) {
        this.parser = parser;
        this.isActiveList = isActiveList;
        this.regionStartLoc = regionStartLoc;
    }

    /**
     * Create a profile of the same class as this object containing just the provided isActiveList
     *
     * Used by clients to create derived activity profiles (such as ones without the starting X
     * sites because they've been removed in an ActiveRegion) of the same class.
     *
     * @param isActiveList the active results list to use in the derived instance
     * @return a freshly allocated data set
     */
    protected ActivityProfile createDerivedProfile(final List<ActivityProfileState> isActiveList) {
        return new ActivityProfile(parser, isActiveList, regionStartLoc);
    }

    @Override
    public String toString() {
        return "ActivityProfile{" +
                "start=" + regionStartLoc +
                ", stop=" + regionStopLoc +
                '}';
    }

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param result a well-formed ActivityProfileState result to incorporate into this profile
     */
    @Requires("result != null")
    public void add(final ActivityProfileState result) {
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

    /**
     * How many profile results are in this profile?
     * @return the number of profile results
     */
    @Ensures("result >= 0")
    public int size() {
        return isActiveList.size();
    }

    /**
     * Is this profile empty?
     * @return true if the profile is empty
     */
    @Ensures("isEmpty() == (size() == 0)")
    public boolean isEmpty() {
        return isActiveList.isEmpty();
    }

    /**
     * Get the list of active profile results in this object
     * @return a non-null, ordered list of active profile results
     */
    @Ensures("result != null")
    protected List<ActivityProfileState> getActiveList() {
        return isActiveList;
    }

    /**
     * Finalize the probabilities in this activity profile, preparing it for a future
     * call to createActiveRegions.  This function returns a new profile with cleaned
     * up activity estimates.
     *
     * This code looks at the current list of states, cleans them up, and then returns
     * a newly allocated ActivityProfile
     *
     * @return a newly allocated ActivityProfile based on the current state of this
     * profile, but that has been "finalized" as required by the profile implementation
     */
    public ActivityProfile finalizeProfile() {
        int iii = 0;
        for( final double prob : finalizeProbabilities() ) {
            final ActivityProfileState result = isActiveList.get(iii++);
            result.isActiveProb = prob;
            result.resultState = ActivityProfileState.Type.NONE;
            result.resultValue = null;
        }

        return createDerivedProfile(isActiveList);
    }

    public double[] finalizeProbabilities() {
        final double[] activeProbArray = new double[isActiveList.size()];

        int iii = 0;
        for( final ActivityProfileState result : isActiveList ) {
            activeProbArray[iii++] = result.isActiveProb;
        }

        iii = 0;
        for( final ActivityProfileState result : isActiveList ) {
            if( result.resultState.equals(ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS) ) { // special code to deal with the problem that high quality soft clipped bases aren't added to pileups
                final int numHQClips = result.resultValue.intValue();
                for( int jjj = Math.max(0, iii - numHQClips); jjj < Math.min(activeProbArray.length, iii+numHQClips); jjj++ ) {
                    activeProbArray[jjj] = Math.max(activeProbArray[jjj], activeProbArray[iii]);
                }
            }
            iii++;
        }

        return activeProbArray;
    }

    /**
     * Partition this profile into active regions
     * @param activeRegionExtension the amount of margin overlap in the active region
     * @return the list of active regions
     */
    public List<ActiveRegion> createActiveRegions( final int activeRegionExtension, final int maxRegionSize ) {
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
    private List<ActiveRegion> createActiveRegion(final boolean isActive, final int curStart, final int curEnd, final int activeRegionExtension, final int maxRegionSize) {
        return createActiveRegion(isActive, curStart, curEnd, activeRegionExtension, maxRegionSize, new ArrayList<ActiveRegion>());
    }

    private List<ActiveRegion> createActiveRegion(final boolean isActive, final int curStart, final int curEnd, final int activeRegionExtension, final int maxRegionSize, final List<ActiveRegion> returnList) {
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
