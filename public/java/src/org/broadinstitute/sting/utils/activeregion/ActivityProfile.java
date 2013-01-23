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

import java.util.*;

/**
 * Class holding information about per-base activity scores for the
 * active region traversal
 *
 * @author Mark DePristo
 * @since Date created
 */
public class ActivityProfile {
    private final static int MAX_PROB_PROPOGATION_DISTANCE = 10;
    private final static double ACTIVE_PROB_THRESHOLD = 0.002; // TODO: needs to be set-able by the walker author

    protected final List<ActivityProfileState> stateList;
    protected final GenomeLocParser parser;

    protected GenomeLoc regionStartLoc = null;
    protected GenomeLoc regionStopLoc = null;

    /**
     * Create a new empty ActivityProfile
     * @param parser the parser we can use to create genome locs, cannot be null
     */
    public ActivityProfile(final GenomeLocParser parser) {
        if ( parser == null ) throw new IllegalArgumentException("parser cannot be null");

        this.parser = parser;
        this.stateList = new ArrayList<ActivityProfileState>();
    }

    @Override
    public String toString() {
        return "ActivityProfile{" +
                "start=" + regionStartLoc +
                ", stop=" + regionStopLoc +
                '}';
    }

    /**
     * How far away can probability mass be moved around in this profile?
     *
     * This distance puts an upper limit on how far, in bp, we will ever propogate probability max around
     * when adding a new ActivityProfileState.  For example, if the value of this function is
     * 10, and you are looking at a state at bp 5, and we know that no states beyond 5 + 10 will have
     * their probability propograted back to that state.
     *
     * @return a positive integer distance in bp
     */
    @Ensures("result >= 0")
    public int getMaxProbPropagationDistance() {
        return MAX_PROB_PROPOGATION_DISTANCE;
    }

    /**
     * How many profile results are in this profile?
     * @return the number of profile results
     */
    @Ensures("result >= 0")
    public int size() {
        return stateList.size();
    }

    /**
     * Is this profile empty?
     * @return true if the profile is empty
     */
    @Ensures("isEmpty() == (size() == 0)")
    public boolean isEmpty() {
        return stateList.isEmpty();
    }

    /**
     * Get the span of this activity profile, which is from the start of the first state to the stop of the last
     * @return a potentially null GenomeLoc.  Will be null if this profile is empty
     */
    public GenomeLoc getSpan() {
        return isEmpty() ? null : regionStartLoc.endpointSpan(regionStopLoc);
    }

    @Requires("! isEmpty()")
    public int getContigIndex() {
        return regionStartLoc.getContigIndex();
    }

    @Requires("! isEmpty()")
    public int getStop() {
        return regionStopLoc.getStop();
    }

    /**
     * Get the list of active profile results in this object
     * @return a non-null, ordered list of active profile results
     */
    @Ensures("result != null")
    protected List<ActivityProfileState> getStateList() {
        return stateList;
    }

    /**
     * Get the probabilities of the states as a single linear array of doubles
     * @return a non-null array
     */
    @Ensures("result != null")
    protected double[] getProbabilitiesAsArray() {
        final double[] probs = new double[getStateList().size()];
        int i = 0;
        for ( final ActivityProfileState state : getStateList() )
            probs[i++] = state.isActiveProb;
        return probs;
    }

    /**
     * Helper function that gets the genome loc for a site offset from relativeLoc, protecting ourselves from
     * falling off the edge of the contig.
     *
     * @param relativeLoc the location offset is relative to
     * @param offset the offset from relativeLoc where we'd like to create a GenomeLoc
     * @return a genome loc with relativeLoc.start + offset, if this is on the contig, null otherwise
     */
    @Requires("relativeLoc != null")
    protected GenomeLoc getLocForOffset(final GenomeLoc relativeLoc, final int offset) {
        final int start = relativeLoc.getStart() + offset;
        if ( start < 0 || start > getCurrentContigLength() ) {
            return null;
        } else {
            return parser.createGenomeLoc(regionStartLoc.getContig(), start);
        }
    }

    /**
     * Get the length of the current contig
     * @return the length in bp
     */
    @Requires("regionStartLoc != null")
    @Ensures("result > 0")
    private int getCurrentContigLength() {
        // TODO -- fix performance problem with getContigInfo
        return parser.getContigInfo(regionStartLoc.getContig()).getSequenceLength();
    }

    // --------------------------------------------------------------------------------
    //
    // routines to add states to a profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    @Requires("state != null")
    public void add(final ActivityProfileState state) {
        final GenomeLoc loc = state.getLoc();

        if ( regionStartLoc == null ) {
            regionStartLoc = loc;
            regionStopLoc = loc;
        } else {
            // TODO -- need to figure out where to add loc as the regions will be popping off the front
            if ( regionStopLoc.getStart() != loc.getStart() - 1 )
                throw new IllegalArgumentException("Bad add call to ActivityProfile: loc " + loc + " not immediate after last loc " + regionStopLoc );
            regionStopLoc = loc;
        }

        final Collection<ActivityProfileState> processedStates = processState(state);
        for ( final ActivityProfileState processedState : processedStates ) {
            incorporateSingleState(processedState);
        }
    }

    /**
     * Incorporate a single activity profile state into the current list of states
     *
     * If state's position occurs immediately after the last position in this profile, then
     * the state is appended to the state list.  If it's within the existing states list,
     * the prob of stateToAdd is added to its corresponding state in the list.  If the
     * position would be before the start of this profile, stateToAdd is simply ignored.
     *
     * @param stateToAdd the state we want to add to the states list
     */
    @Requires("stateToAdd != null")
    private void incorporateSingleState(final ActivityProfileState stateToAdd) {
        final int position = stateToAdd.getOffset(regionStartLoc);

        if ( position > size() )
            // should we allow this?  probably not
            throw new IllegalArgumentException("Must add state contiguous to existing states");

        if ( position >= 0 ) {
            // ignore states starting before this regions start
            if ( position < size() ) {
                stateList.get(position).isActiveProb += stateToAdd.isActiveProb;
            } else {
                if ( position != size() ) throw new IllegalStateException("position == size but it wasn't");
                stateList.add(stateToAdd);
            }
        }
    }

    /**
     * Process justAddedState, returning a collection of derived states that actually be added to the stateList
     *
     * The purpose of this function is to transform justAddedStates, if needed, into a series of atomic states
     * that we actually want to track.  For example, if state is for soft clips, we transform that single
     * state into a list of states that surround the state up to the distance of the soft clip.
     *
     * Can be overridden by subclasses to transform states in any way
     *
     * There's no particular contract for the output states, except that they can never refer to states
     * beyond the current end of the stateList unless the explictly include preceding states before
     * the reference.  So for example if the current state list is [1, 2, 3] this function could return
     * [1,2,3,4,5] but not [1,2,3,5].
     *
     * @param justAddedState the state our client provided to use to add to the list
     * @return a list of derived states that should actually be added to this profile's state list
     */
    protected Collection<ActivityProfileState> processState(final ActivityProfileState justAddedState) {
        if ( justAddedState.resultState.equals(ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS) ) {
            // special code to deal with the problem that high quality soft clipped bases aren't added to pileups
            final List<ActivityProfileState> states = new LinkedList<ActivityProfileState>();
            final int numHQClips = justAddedState.resultValue.intValue();
            for( int jjj = - numHQClips; jjj <= numHQClips; jjj++ ) {
                final GenomeLoc loc = getLocForOffset(justAddedState.getLoc(), jjj);
                if ( loc != null )
                    states.add(new ActivityProfileState(loc, justAddedState.isActiveProb));
            }

            return states;
        } else {
            return Collections.singletonList(justAddedState);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // routines to get active regions from the profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Get the next completed active regions from this profile, and remove all states supporting them from this profile
     *
     * Takes the current profile and finds all of the active / inactive from the start of the profile that are
     * ready.  By ready we mean unable to have their probability modified any longer by future additions to the
     * profile.  The regions that are popped off the profile take their states with them, so the start of this
     * profile will always be after the end of the last region returned here.
     *
     * The regions are returned sorted by genomic position.
     *
     * This function may not return anything in the list, if no regions are ready
     *
     * No returned region will be larger than maxRegionSize.
     *
     * @param activeRegionExtension the extension value to provide to the constructed regions
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return a non-null list of active regions
     */
    @Ensures("result != null")
    public List<ActiveRegion> popReadyActiveRegions(final int activeRegionExtension, final int maxRegionSize, final boolean forceConversion) {
        if ( activeRegionExtension < 0 ) throw new IllegalArgumentException("activeRegionExtension must be >= 0 but got " + activeRegionExtension);
        if ( maxRegionSize < 1 ) throw new IllegalArgumentException("maxRegionSize must be >= 1 but got " + maxRegionSize);

        final LinkedList<ActiveRegion> regions = new LinkedList<ActiveRegion>();

        while ( true ) {
            final ActiveRegion nextRegion = popNextReadyActiveRegion(activeRegionExtension, maxRegionSize, forceConversion);
            if ( nextRegion == null )
                return regions;
            else {
                regions.add(nextRegion);
            }
        }
    }

    /**
     * Helper function for popReadyActiveRegions that pops the first ready region off the front of this profile
     *
     * If a region is returned, modifies the state of this profile so that states used to make the region are
     * no longer part of the profile.  Associated information (like the region start position) of this profile
     * are also updated.
     *
     * @param activeRegionExtension the extension value to provide to the constructed regions
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return a fully formed active region, or null if none can be made
     */
    private ActiveRegion popNextReadyActiveRegion(final int activeRegionExtension, final int maxRegionSize, final boolean forceConversion) {
        if ( stateList.isEmpty() )
            return null;

        final ActivityProfileState first = stateList.get(0);
        final boolean isActiveRegion = first.isActiveProb > ACTIVE_PROB_THRESHOLD;
        final int offsetOfNextRegionEnd = findEndOfRegion(isActiveRegion, maxRegionSize, forceConversion);
        if ( offsetOfNextRegionEnd == -1 )
            // couldn't find a valid ending offset, so we return null
            return null;

        // we need to create the active region, and clip out the states we're extracting from this profile
        final List<ActivityProfileState> sub = stateList.subList(0, offsetOfNextRegionEnd + 1);
        final List<ActivityProfileState> supportingStates = new ArrayList<ActivityProfileState>(sub);
        sub.clear();

        // update the start and stop locations as necessary
        if ( stateList.isEmpty() ) {
            regionStartLoc = regionStopLoc = null;
        } else {
            regionStartLoc = stateList.get(0).getLoc();
        }
        final GenomeLoc regionLoc = parser.createGenomeLoc(first.getLoc().getContig(), first.getLoc().getStart(), first.getLoc().getStart() + offsetOfNextRegionEnd);
        return new ActiveRegion(regionLoc, supportingStates, isActiveRegion, parser, activeRegionExtension);
    }

    /**
     * Find the end of the current region, returning the index into the element isActive element, or -1 if the region isn't done
     *
     * The current region is defined from the start of the stateList, looking for elements that have the same isActiveRegion
     * flag (i.e., if isActiveRegion is true we are looking for states with isActiveProb > threshold, or alternatively
     * for states < threshold).  The maximize size of the returned region is maxRegionSize.  If forceConversion is
     * true, then we'll return the region end even if this isn't safely beyond the max prob propogation distance.
     *
     * @param isActiveRegion is the region we're looking for an active region or inactive region?
     * @param maxRegionSize the maximize size of the returned region
     * @param forceConversion if true, we'll return a region whose end isn't sufficiently far from the end of the
     *                        stateList.  Used to close out the active region when we've hit some kind of end (such
     *                        as the end of the contig)
     * @return the index into stateList of the last element of this region, or -1 if it cannot be found
     */
    @Ensures({
            "result >= -1",
            "result == -1 || result < maxRegionSize",
            "! (result == -1 && forceConversion)"})
    private int findEndOfRegion(final boolean isActiveRegion, final int maxRegionSize, final boolean forceConversion) {
        int i = 0;
        while ( i < stateList.size() && i < maxRegionSize ) {
            if ( stateList.get(i).isActiveProb > ACTIVE_PROB_THRESHOLD != isActiveRegion ) {
                break;
            }
            i++;
        }

        // we're one past the end, so i must be decremented
        return forceConversion || i + getMaxProbPropagationDistance() < stateList.size() ? i - 1 : -1;
    }
}
