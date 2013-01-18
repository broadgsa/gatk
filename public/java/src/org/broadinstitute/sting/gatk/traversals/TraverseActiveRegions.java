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

package org.broadinstitute.sting.gatk.traversals;

import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.progressmeter.ProgressMeter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Implement active region traversal
 *
 * User: depristo
 * Date: 1/9/13
 * Time: 4:45 PM
 *
 * Live region:
 *
 *   The ART tracks a thing called the live region.  The live region is a position on a specific contig
 *   of the alignment start of the last read we processed during this traversal.  Because the
 *   read stream is sorted, future reads must occurs in the the live region.  Therefore the the dead region
 *   (everything to the left of the live boundary) cannot have any more read data.  The live / dead
 *   regions are used to decide when we can safely call map on active regions, as only active regions
 *   contained completely within the dead region (including extensions) have a complete set of read data
 *   in the collected read list.  All of the data related to the live region is captured by the local
 *   variable spanOfLastReadSeen
 *
 */
public class TraverseActiveRegions<M, T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,LocusShardDataProvider> {
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);
    protected final static boolean DEBUG = false;

    // set by the tranversal
    private int activeRegionExtension = -1;
    private int maxRegionSize = -1;

    private final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();

    private LinkedList<GATKSAMRecord> myReads = new LinkedList<GATKSAMRecord>();
    private GenomeLoc spanOfLastReadSeen = null;

    protected int getActiveRegionExtension() {
        return activeRegionExtension;
    }

    protected int getMaxRegionSize() {
        return maxRegionSize;
    }

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public String toString() {
        return "TraverseActiveRegions";
    }

    @Override
    public void initialize(GenomeAnalysisEngine engine, Walker walker, ProgressMeter progressMeter) {
        super.initialize(engine, walker, progressMeter);
        activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

        final ActiveRegionWalker arWalker = (ActiveRegionWalker)walker;
        if ( arWalker.wantsExtendedReads() && ! arWalker.wantsNonPrimaryReads() ) {
            throw new IllegalArgumentException("Active region walker " + arWalker + " requested extended events but not " +
                    "non-primary reads, an inconsistent state.  Please modify the walker");
        }
    }

    /**
     * Is the loc outside of the intervals being requested for processing by the GATK?
     * @param loc
     * @return
     */
    protected boolean outsideEngineIntervals(final GenomeLoc loc) {
        return engine.getIntervals() != null && ! engine.getIntervals().overlaps(loc);
    }

    /**
     * Take the individual isActive calls and integrate them into contiguous active regions and
     * add these blocks of work to the work queue
     * band-pass filter the list of isActive probabilities and turn into active regions
     *
     * @param profile
     * @param activeRegions
     * @return
     */
    protected ActivityProfile incorporateActiveRegions(final ActivityProfile profile,
                                                       final List<ActiveRegion> activeRegions) {
        if ( profile.isEmpty() )
            throw new IllegalStateException("trying to incorporate an empty active profile " + profile);

        final ActivityProfile bandPassFiltered = profile.bandPassFilter();
        activeRegions.addAll(bandPassFiltered.createActiveRegions( getActiveRegionExtension(), getMaxRegionSize() ));
        return new ActivityProfile( engine.getGenomeLocParser(), profile.hasPresetRegions() );
    }

    protected final ActivityProfileResult walkerActiveProb(final ActiveRegionWalker<M, T> walker,
                                                           final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                           final AlignmentContext locus, final GenomeLoc location) {
        if ( walker.hasPresetActiveRegions() ) {
            return new ActivityProfileResult(location, walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0);
        } else {
            return walker.isActive( tracker, refContext, locus );
        }
    }

    protected ReferenceOrderedView getReferenceOrderedView(final ActiveRegionWalker<M, T> walker,
                                                           final LocusShardDataProvider dataProvider,
                                                           final LocusView locusView) {
        if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
            return new ManagingReferenceOrderedView( dataProvider );
        else
            return (RodLocusView)locusView;
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param walker
     */
    protected void writeActiveRegionsToStream( final ActiveRegionWalker<M, T> walker ) {
        // Just want to output the active regions to a file, not actually process them
        for( final ActiveRegion activeRegion : workQueue ) {
            if( activeRegion.isActive ) {
                walker.activeRegionOutStream.println( activeRegion.getLocation() );
            }
        }
    }

    /**
     * Did read appear in the last shard?
     *
     * When we transition across shard boundaries we see duplicate reads because
     * each shard contains the reads that *overlap* the shard.  So if we just finished
     * shard 1-1000 and are now in 1001-2000 we'll see duplicate reads from 1001
     * that overlapped 1-1000.  This function tests read to determine if we would have
     * seen it before by asking if read.getAlignmentStart() is less than the
     * stop position of the last seen read at the start of the traversal.  The reason
     * we need to use the location of the last read at the start of the traversal
     * is that we update the lastRead during the traversal, and we only want to filter
     * out reads whose start is before the last read of the previous shard, not the
     * current shard.
     *
     * @param locOfLastReadAtTraversalStart the location of the last read seen at the start of the traversal
     * @param read the read we want to test if it's already been seen in the last shard
     * @return true if read would have appeared in the last shard, false otherwise
     */
    protected boolean appearedInLastShard(final GenomeLoc locOfLastReadAtTraversalStart, final GATKSAMRecord read) {
        if ( locOfLastReadAtTraversalStart == null )
            // we're in the first shard, so obviously the answer is no
            return false;
        else {
            // otherwise check to see if the alignment occurred in the previous shard
            return read.getAlignmentStart() <= locOfLastReadAtTraversalStart.getStart()
                    // we're on the same contig
                    && read.getReferenceIndex() == locOfLastReadAtTraversalStart.getContigIndex();
        }

    }

    // -------------------------------------------------------------------------------------
    //
    // Actual traverse function
    //
    // -------------------------------------------------------------------------------------

    /**
     * Is the current shard on a new contig w.r.t. the previous shard?
     * @param currentShard the current shard we are processing
     * @return true if the last shard was on a different contig than the current shard
     */
    private boolean onNewContig(final Shard currentShard) {
        return spanOfLastSeenRead() != null
                && spanOfLastSeenRead().getContigIndex() != currentShard.getLocation().getContigIndex();
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseActiveRegions.traverse: Shard is %s", dataProvider));

        final LocusView locusView = new AllLocusView(dataProvider);

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

        final List<ActiveRegion> activeRegions = new LinkedList<ActiveRegion>();
        ActivityProfile profile = new ActivityProfile(engine.getGenomeLocParser(), walker.hasPresetActiveRegions() );

        ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, dataProvider, locusView);

        // We keep processing while the next reference location is within the interval
        final GenomeLoc locOfLastReadAtTraversalStart = spanOfLastSeenRead();

        // if we've moved onto a new contig, process all of the active regions
        if ( onNewContig(dataProvider.getShard()) )
            sum = processActiveRegions(walker, sum, true);

        GenomeLoc prevLoc = null;
        while( locusView.hasNext() ) {
            final AlignmentContext locus = locusView.next();
            final GenomeLoc location = locus.getLocation();

            // Grab all the previously unseen reads from this pileup and add them to the massive read list
            // Note that this must occur before we leave because we are outside the intervals because
            // reads may occur outside our intervals but overlap them in the future
            final Collection<GATKSAMRecord> reads = locusView.getLIBS().transferReadsFromAllPreviousPileups();
            for( final GATKSAMRecord read : reads ) {
                if ( appearedInLastShard(locOfLastReadAtTraversalStart, read) ) {
                    if ( DEBUG ) logger.warn("Skipping duplicated " + read.getReadName());
                } else {
                    if ( DEBUG ) logger.warn("Adding read " + read.getReadName() + " at " + engine.getGenomeLocParser().createGenomeLoc(read) + " from provider " + dataProvider);
                    rememberLastReadLocation(read);
                    myReads.add(read);
                }
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

        // now go and process all of the active regions
        sum = processActiveRegions(walker, sum, false);

        return sum;
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal(final Walker<M, T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M, T>)walker, sum, true);
    }

    // -------------------------------------------------------------------------------------
    //
    // Functions to manage and interact with the live / dead zone
    //
    // -------------------------------------------------------------------------------------

    /**
     * Update the live region to reflect that the last read we've seen in the traversal is read
     *
     * Requires that sequential calls always be provided reads in coordinate sorted order
     *
     * @param read the last read we've seen during the traversal
     */
    protected void rememberLastReadLocation(final GATKSAMRecord read) {
        final GenomeLoc currentLocation = engine.getGenomeLocParser().createGenomeLoc(read);
        if ( spanOfLastReadSeen == null )
            spanOfLastReadSeen = currentLocation;
        else {
            if ( currentLocation.isBefore(spanOfLastReadSeen) )
                throw new IllegalStateException("Updating last read seen in the traversal with read " + read + " with span " + currentLocation + " but this occurs before the previously seen read " + spanOfLastReadSeen);
            spanOfLastReadSeen = currentLocation;
        }
    }

    /**
     * Get a GenomeLoc indicating the start (heading to the right) of the live ART region.
     * @return the left-most position of the live region on the genome
     */
    protected GenomeLoc spanOfLastSeenRead() {
        return spanOfLastReadSeen;
    }

    /**
     * Is the active region completely within the traversal's dead zone?
     *
     * @param region the region we want to test
     * @return true if the extended location of region is completely within the current dead zone, false otherwise
     */
    protected boolean regionCompletelyWithinDeadZone(final ActiveRegion region) {
        return region.getExtendedLoc().getStop() < spanOfLastSeenRead().getStart()
                || ! region.getExtendedLoc().onSameContig(spanOfLastSeenRead());
    }

    /**
     * Is the read dead?  That is, can it no longer be in any future active region, and therefore can be discarded?
     *
     * read: start |--------> stop ------ stop + extension
     * region:                      start |-----------------| end
     *
     * Since the regions are coming in order, read could potentially be contained in a future interval if
     * stop + activeRegionExtension >= end.  If, on the other hand, stop + extension is < the end
     * of this region, then we can discard it, since any future region could only include reads
     * up to end + 1 - extension.
     *
     * Note that this function doesn't care about the dead zone.  We're assuming that by
     * actually calling this function with an active region that region is already in the dead zone,
     * so checking that the read is in the dead zone doesn't make sense.
     *
     * @param read the read we're testing
     * @param activeRegion the current active region
     * @return true if the read is dead, false other
     */
    @Requires({"read != null", "activeRegion != null"})
    private boolean readCannotOccurInAnyMoreActiveRegions(final GATKSAMRecord read, final ActiveRegion activeRegion) {
        return read.getAlignmentEnd() + getActiveRegionExtension() < activeRegion.getLocation().getStop();
    }

    // -------------------------------------------------------------------------------------
    //
    // Functions to process active regions that are ready for map / reduce calls
    //
    // -------------------------------------------------------------------------------------

    private T processActiveRegions(final ActiveRegionWalker<M, T> walker, T sum, final boolean forceRegionsToBeActive) {
        if( walker.activeRegionOutStream != null ) {
            writeActiveRegionsToStream(walker);
            return sum;
        } else {
            return callWalkerMapOnActiveRegions(walker, sum, forceRegionsToBeActive);
        }
    }

    private T callWalkerMapOnActiveRegions(final ActiveRegionWalker<M, T> walker, T sum, final boolean forceRegionsToBeActive) {
        // Since we've traversed sufficiently past this point (or this contig!) in the workQueue we can unload those regions and process them
        // TODO can implement parallel traversal here
        while( workQueue.peek() != null ) {
            final ActiveRegion activeRegion = workQueue.peek();
            if ( forceRegionsToBeActive || regionCompletelyWithinDeadZone(activeRegion) ) {
                if ( DEBUG ) logger.warn("Processing active region " + activeRegion + " dead zone " + spanOfLastSeenRead());
                sum = processActiveRegion( workQueue.remove(), sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    protected T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker) {
        final Iterator<GATKSAMRecord> liveReads = myReads.iterator();
        while ( liveReads.hasNext() ) {
            boolean killed = false;
            final GATKSAMRecord read = liveReads.next();
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );

            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                activeRegion.add(read);

                if ( ! walker.wantsNonPrimaryReads() ) {
                    if ( DEBUG ) logger.warn("Removing read " + read.getReadName() + " at " + readLoc + " with dead zone start " + spanOfLastSeenRead());
                    liveReads.remove();
                    killed = true;
                }
            } else if( walker.wantsExtendedReads() && activeRegion.getExtendedLoc().overlapsP( readLoc )) {
                activeRegion.add( read );
            }

            if ( ! killed && readCannotOccurInAnyMoreActiveRegions(read, activeRegion) ) {
                if ( DEBUG ) logger.warn("Removing read " + read.getReadName() + " at " + readLoc + " with dead zone start " + spanOfLastSeenRead());
                liveReads.remove();
            }
        }

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLoc());
        final M x = walker.map(activeRegion, null);
        return walker.reduce( x, sum );
    }
}
