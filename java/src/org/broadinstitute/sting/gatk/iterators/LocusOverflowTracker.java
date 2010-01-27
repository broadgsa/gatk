package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

/**
 * a helper class that organizes the output of warning messages from read pile-ups that
 * are greater than the max pile-up size.  We only want a single warning from each non-contigous
 * interval, up until the maximum warning limit.
 *
 * cleanWarningQueue() must be called when you're finished with the LocusOverflowTracker to make
 * sure that no errors are left in the queue.
 *
 */
public class LocusOverflowTracker {
    // the last interval we emitted a warning for
    protected GenomeLoc lastLocation = null;

    // the maximum warning count, and the number of warnings emitted
    protected static int warningsEmitted = 0;
    public static final int MAX_WARNINGS = 100;

    // our maximum pileup size
    protected final int maxPileupSize;

    // do we have a pending warning?
    protected boolean warningInQueue = false;

    /**
     * create a LocusOverflowTracker
     *
     * @param maxPileup the maximum allowed pile-up size
     */
    public LocusOverflowTracker(int maxPileup) {
        warningInQueue = false;
        maxPileupSize = maxPileup;
    }

    /**
     * have we exceeded the maximum pile-up size?
     *
     * @param loc        the current location
     * @param pileupSize the pile-up size
     *
     * @return return true if we're greater, false if we're not
     */
    public boolean exceeded(GenomeLoc loc, int pileupSize) {
        boolean exceeded = pileupSize >= maxPileupSize;
        if (exceeded) {

            // if the last location is null, we're starting a new region that exceeds max_reads_at_locus
            if (lastLocation == null) lastLocation = loc;
            // are we contiguous to the last genome loc?
            else if (lastLocation.contiguousP(loc)) {
                lastLocation = lastLocation.merge(loc);
            }
            // we have an existing region, and the current is not contiguous.  Emit the old and store the new
            else {
                warnUser();
                lastLocation = loc;
            }

            // regardless, we have a warning in the queue
            warningInQueue = true;
        }
        // we don't have a warning at this, but there is one in the queue
        else if (warningInQueue) {
            warnUser();
            lastLocation = null;
        }
        // return true if we exceeded the max size at this location
        return exceeded;
    }

    /**
     * clean up the warning queue, making sure we haven't stored a warning
     * that hasn't been emitted yet.
     */
    public void cleanWarningQueue() {
        if (warningInQueue) warnUser();
    }

    /** warn the user, checking to make sure we haven't exceeded the maximum warning level. */
    protected void warnUser() {

        // make sure we have a warning in the queue
        if (!warningInQueue) throw new IllegalStateException("Without a warning in the queue, we shouldn't see a call to warnUser()");

        // reset the warning light
        warningInQueue = false;

        // check to see if we've meet our warning threshold or not. If we're equal to the threshold emit a message saying this
        // is the last warning they'll see
        if (warningsEmitted < MAX_WARNINGS) {
            warningsEmitted++;
            Utils.warnUser("Unable to add reads to the pile-up, we're over the hanger limit of " + maxPileupSize + " at location: " + lastLocation);
        } else if (warningsEmitted == MAX_WARNINGS) {
            warningsEmitted++;
            Utils.warnUser("Unable to add reads to the pile-up, we're over the hanger limit of " + maxPileupSize + " at location: " + lastLocation +
                    "; the maximum warning count has been reached, we will no longer emit warnings of this nature!!");
        }
    }

    /**
     * is the specified location in the current exceeded pileup region
     * @param position the position
     * @return true if we're in that region
     */
    public boolean inDroppedRegion(GenomeLoc position) {
        if (lastLocation == null || position == null) return false;
        return position.overlapsP(lastLocation) ? true : false;
    }
}
