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

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.ArrayList;
import java.util.Arrays;
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
    final List<Double> isActiveList;

    // todo -- add upfront the start and stop of the intervals
    // todo -- check that no regions are unexpectedly missing
    // todo -- add unit tests
    // TODO -- own preset regions
    public ActivityProfile(final GenomeLocParser parser, final boolean presetRegions) {
        this(parser, presetRegions, new ArrayList<Double>(), null);
    }

    protected ActivityProfile(final GenomeLocParser parser, final boolean presetRegions, final List<Double> isActiveList, final GenomeLoc regionStartLoc) {
        this.parser = parser;
        this.presetRegions = presetRegions;
        this.isActiveList = isActiveList;
        this.regionStartLoc = regionStartLoc;
    }

    public void add(final GenomeLoc loc, final double score) {
        // todo -- test for validity
        isActiveList.add(score);
        if( regionStartLoc == null ) {
            regionStartLoc = loc;
        }
    }

    public int size() {
        return isActiveList.size();
    }

    /**
     * Band pass this ActivityProfile, producing a new profile that's band pass filtered
     * @return a new ActivityProfile that's the band-pass filtered version of this profile
     */
    public ActivityProfile bandPassFilter() {
        final Double[] activeProbArray = isActiveList.toArray(new Double[isActiveList.size()]);
        final Double[] filteredProbArray = new Double[activeProbArray.length];
        final int FILTER_SIZE = ( presetRegions ? 0 : 50 ); // TODO: needs to be set-able by the walker author
        for( int iii = 0; iii < activeProbArray.length; iii++ ) {
            double maxVal = 0;
            for( int jjj = Math.max(0, iii-FILTER_SIZE); jjj < Math.min(isActiveList.size(), iii+FILTER_SIZE+1); jjj++ ) {
                if( activeProbArray[jjj] > maxVal ) { maxVal = activeProbArray[jjj]; }
            }
            filteredProbArray[iii] = maxVal;
        }

        return new ActivityProfile(parser, presetRegions, Arrays.asList(filteredProbArray), regionStartLoc);
    }

    /**
     * Partition this profile into active regions
     * @param activeRegionExtension
     * @return
     */
    public List<ActiveRegion> createActiveRegions( final int activeRegionExtension ) {
        final int MAX_ACTIVE_REGION = ( presetRegions ? 16001 : 425 ); // TODO: needs to be set-able by the walker author
        final double ACTIVE_PROB_THRESHOLD = 0.2; // TODO: needs to be set-able by the walker author

        if( isActiveList.size() == 0 ) {
            // no elements in the active list, just return an empty one
            return Collections.emptyList();
        } else if( isActiveList.size() == 1 ) {
            // there's a single element, it's either active or inactive
            boolean isActive = isActiveList.get(0) > ACTIVE_PROB_THRESHOLD;
            final ActiveRegion region = createActiveRegion(isActive, 0, 0, activeRegionExtension );
            return Collections.singletonList(region);
        } else {
            // there are 2+ elements, divide these up into regions
            final ArrayList<ActiveRegion> returnList = new ArrayList<ActiveRegion>();
            boolean isActive = isActiveList.get(0) > ACTIVE_PROB_THRESHOLD;
            int curStart = 0;
            for(int iii = 1; iii < isActiveList.size(); iii++ ) {
                final boolean thisStatus = isActiveList.get(iii) > ACTIVE_PROB_THRESHOLD;
                if( isActive != thisStatus || (iii-curStart) > MAX_ACTIVE_REGION ) {
                    returnList.add( createActiveRegion(isActive, curStart, iii-1, activeRegionExtension) );
                    isActive = thisStatus;
                    curStart = iii;
                }
            }

            if( curStart != isActiveList.size()-1 ) {
                returnList.add( createActiveRegion(isActive, curStart, isActiveList.size()-1, activeRegionExtension) );
            }
            return returnList;
        }
    }

    /**
     * Helper routine to create an active region based on our current start and end offsets
     * @param isActive should the region be active?
     * @param curStart offset (0-based) from the start of this region
     * @param curEnd offset (0-based) from the start of this region
     * @param activeRegionExtension
     * @return a fully initialized ActiveRegion with the above properties
     */
    private final ActiveRegion createActiveRegion(final boolean isActive, final int curStart, final int curEnd, final int activeRegionExtension) {
        final GenomeLoc loc = parser.createGenomeLoc(regionStartLoc.getContig(), regionStartLoc.getStart() + curStart, regionStartLoc.getStart() + curEnd);
        return new ActiveRegion( loc, isActive, parser, activeRegionExtension );
    }
}
