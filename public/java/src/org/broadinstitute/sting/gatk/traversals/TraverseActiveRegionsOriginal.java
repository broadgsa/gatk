package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/9/11
 */

public class TraverseActiveRegionsOriginal<M,T> extends TraverseActiveRegions<M,T> {
    private final LinkedHashSet<GATKSAMRecord> myReads = new LinkedHashSet<GATKSAMRecord>();

    protected Collection<GATKSAMRecord> getReadsInCurrentRegion() {
        return myReads;
    }

    protected void removeReadsFromCurrentRegion(final List<GATKSAMRecord> placedReads) {
        myReads.removeAll( placedReads ); // remove all the reads which have been placed into their active region
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseActiveRegions.traverse: Shard is %s", dataProvider));

        final LocusView locusView = new AllLocusView(dataProvider);

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

        int minStart = Integer.MAX_VALUE;
        final List<ActiveRegion> activeRegions = new LinkedList<ActiveRegion>();
        ActivityProfile profile = new ActivityProfile(engine.getGenomeLocParser(), walker.hasPresetActiveRegions() );

        ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, dataProvider, locusView);

        // We keep processing while the next reference location is within the interval
        GenomeLoc prevLoc = null;
        while( locusView.hasNext() ) {
            final AlignmentContext locus = locusView.next();
            final GenomeLoc location = locus.getLocation();

            // Grab all the previously unseen reads from this pileup and add them to the massive read list
            // Note that this must occur before we leave because we are outside the intervals because
            // reads may occur outside our intervals but overlap them in the future
            // TODO -- this whole HashSet logic should be changed to a linked list of reads with
            // TODO -- subsequent pass over them to find the ones overlapping the active regions
            for( final PileupElement p : locus.getBasePileup() ) {
                final GATKSAMRecord read = p.getRead();
                if( !myReads.contains(read) ) {
                    myReads.add(read);
                }

                // If this is the last pileup for this shard calculate the minimum alignment start so that we know
                // which active regions in the work queue are now safe to process
                minStart = Math.min(minStart, read.getAlignmentStart());
            }

            // skip this location -- it's not part of our engine intervals
            if ( outsideEngineIntervals(location) )
                continue;

            if ( prevLoc != null && location.getStart() != prevLoc.getStop() + 1 ) {
                // we've move across some interval boundary, restart profile
                profile = incorporateActiveRegions(profile, activeRegions);
            }

            dataProvider.getShard().getReadMetrics().incrementNumIterations();

            // create reference context. Note that if we have a pileup of "extended events", the context will
            // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
            final ReferenceContext refContext = referenceView.getReferenceContext(location);

            // Iterate forward to get all reference ordered data covering this location
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

            // Call the walkers isActive function for this locus and add them to the list to be integrated later
            profile.add(walkerActiveProb(walker, tracker, refContext, locus, location));

            prevLoc = location;

            printProgress(locus.getLocation());
        }

        updateCumulativeMetrics(dataProvider.getShard());

        if ( ! profile.isEmpty() )
            incorporateActiveRegions(profile, activeRegions);

        // add active regions to queue of regions to process
        // first check if can merge active regions over shard boundaries
        if( !activeRegions.isEmpty() ) {
            if( !workQueue.isEmpty() ) {
                final ActiveRegion last = workQueue.getLast();
                final ActiveRegion first = activeRegions.get(0);
                if( last.isActive == first.isActive && last.getLocation().contiguousP(first.getLocation()) && last.getLocation().size() + first.getLocation().size() <= getMaxRegionSize() ) {
                    workQueue.removeLast();
                    activeRegions.remove(first);
                    workQueue.add( new ActiveRegion(last.getLocation().union(first.getLocation()), first.isActive, this.engine.getGenomeLocParser(), getActiveRegionExtension()) );
                }
            }
            workQueue.addAll( activeRegions );
        }

        logger.debug("Integrated " + profile.size() + " isActive calls into " + activeRegions.size() + " regions." );

        // set the dead zone to the min.  This is incorrect but necessary because of the way we handle things in processActiveRegion
        notifyOfCurrentPosition(engine.getGenomeLocParser().createGenomeLoc(dataProvider.getLocus().getContig(), minStart));
        // now go and process all of the active regions
        sum = processActiveRegions(walker, sum, false);

        return sum;
    }

    @Override
    public String toString() {
        return "TraverseActiveRegionsOriginal";
    }

    @Override
    protected T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker) {
        final ArrayList<GATKSAMRecord> placedReads = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord read : getReadsInCurrentRegion() ) {
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );

            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                // The region which the highest amount of overlap is chosen as the primary region for the read (tie breaking is done as right most region)
                final ActiveRegion bestRegion = getBestRegion(activeRegion, readLoc);
                bestRegion.add( read );

                // The read is also added to all other regions in which it overlaps but marked as non-primary

                if( walker.wantsNonPrimaryReads() ) {
                    if( !bestRegion.equals(activeRegion) ) {
                        activeRegion.add( read );
                    }
                    for( final ActiveRegion otherRegionToTest : workQueue ) {
                        if( !bestRegion.equals(otherRegionToTest) ) {
                            // check for non-primary vs. extended
                            if ( otherRegionToTest.getLocation().overlapsP( readLoc ) ) {
                                otherRegionToTest.add( read );
                            } else if ( walker.wantsExtendedReads() && otherRegionToTest.getExtendedLoc().overlapsP( readLoc ) ) {
                                otherRegionToTest.add( read );
                            }
                        }
                    }
                }
                placedReads.add( read );
                // check for non-primary vs. extended
            } else if( walker.wantsExtendedReads() && activeRegion.getExtendedLoc().overlapsP( readLoc )) {
                activeRegion.add( read );
            }
        }

        removeReadsFromCurrentRegion(placedReads);
        // WARNING: This hashset relies on reads being exactly equal when they are placed in the list as when they are removed. So the ActiveRegionWalker can't modify the reads in any way.

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLoc());
        final M x = walker.map(activeRegion, null);
        return walker.reduce( x, sum );
    }
}
