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

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/9/11
 */

public class TraverseActiveRegionsOptimized<M,T> extends TraverseActiveRegions<M,T> {
    private LinkedList<GATKSAMRecord> myReads = new LinkedList<GATKSAMRecord>();
    private Shard lastShard = null;

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        if ( DEBUG ) logger.warn(String.format("TraverseActiveRegions.traverse: Shard is %s", dataProvider));

        final HashSet<GATKSAMRecord> maybeDuplicatedReads = new HashSet<GATKSAMRecord>();
        // TODO -- there's got to be a better way to know this
        if ( lastShard != dataProvider.getShard() ) {
            maybeDuplicatedReads.addAll(myReads);
            logger.info("Crossing shard boundary requires us to check for duplicates against " + maybeDuplicatedReads.size() +  " reads");
            if ( DEBUG ) logger.warn("Clearing myReads");
        }
        lastShard = dataProvider.getShard();

        final LocusView locusView = new AllLocusView(dataProvider);

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

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
            final Collection<SAMRecord> reads = locusView.getLIBS().transferReadsFromAllPreviousPileups();
            for( final SAMRecord read : reads ) {
                notifyOfCurrentPosition((GATKSAMRecord)read);
                // most of the time maybeDuplicatedReads is empty
                // TODO -- I believe that because of the ordering of reads that as soon as we don't find a read in the
                // TODO -- potential list of duplicates we can clear the hashset
                if ( ! maybeDuplicatedReads.isEmpty() && maybeDuplicatedReads.contains(read) ) {
                    if ( DEBUG ) logger.warn("Skipping duplicated " + read.getReadName());
                } else {
                    if ( DEBUG ) logger.warn("Adding read " + read.getReadName() + " at " + engine.getGenomeLocParser().createGenomeLoc(read) + " from provider " + dataProvider);
                    myReads.add((GATKSAMRecord)read);
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

    private GenomeLoc startOfLiveRegion = null;

    protected void notifyOfCurrentPosition(final GATKSAMRecord read) {
        notifyOfCurrentPosition(engine.getGenomeLocParser().createGenomeLoc(read));
    }

    protected void notifyOfCurrentPosition(final GenomeLoc currentLocation) {
        if ( startOfLiveRegion == null )
            startOfLiveRegion = currentLocation;
        else
            startOfLiveRegion = startOfLiveRegion.max(currentLocation.getStartLocation());
    }

    protected GenomeLoc getStartOfLiveRegion() {
        return startOfLiveRegion;
    }

    protected boolean regionCompletelyWithinDeadZone(final GenomeLoc region, final boolean includeExtension) {
        return (region.getStop() < (getStartOfLiveRegion().getStart() - (includeExtension ? getActiveRegionExtension() : 0)))
                || ! region.onSameContig(getStartOfLiveRegion());
    }

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
            final GenomeLoc extendedLoc = workQueue.peek().getExtendedLoc();
            if ( forceRegionsToBeActive || regionCompletelyWithinDeadZone(extendedLoc, false) ) {
                final ActiveRegion activeRegion = workQueue.remove();
                if ( DEBUG ) logger.warn("Processing active region " + activeRegion + " dead zone " + getStartOfLiveRegion());
                sum = processActiveRegion( activeRegion, sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    @Override
    public String toString() {
        return "TraverseActiveRegionsOptimized";
    }

    private boolean readIsDead(final GATKSAMRecord read, final GenomeLoc readLoc, final ActiveRegion activeRegion) {
        return readLoc.getStop() < activeRegion.getLocation().getStart() && regionCompletelyWithinDeadZone(readLoc, true);
    }

    @Override
    protected T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker) {
        final Iterator<GATKSAMRecord> liveReads = myReads.iterator();
        while ( liveReads.hasNext() ) {
            boolean killed = false;
            final GATKSAMRecord read = liveReads.next();
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );

            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                activeRegion.add(read);

                if ( ! walker.wantsNonPrimaryReads() ) {
                    if ( DEBUG ) logger.warn("Removing read " + read.getReadName() + " at " + readLoc + " with dead zone start " + getStartOfLiveRegion());
                    liveReads.remove();
                    killed = true;
                }
            } else if( walker.wantsExtendedReads() && activeRegion.getExtendedLoc().overlapsP( readLoc )) {
                activeRegion.add( read );
            }

            if ( ! killed && readIsDead(read, readLoc, activeRegion) ) {
                if ( DEBUG ) logger.warn("Removing read " + read.getReadName() + " at " + readLoc + " with dead zone start " + getStartOfLiveRegion());
                liveReads.remove();
            }
        }

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLoc());
        final M x = walker.map(activeRegion, null);
        return walker.reduce( x, sum );
    }


    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    @Override
    public T endTraversal(final Walker<M, T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M, T>)walker, sum, true);
    }

}
