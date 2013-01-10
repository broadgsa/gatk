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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 1/9/13
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TraverseActiveRegions<M, T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,LocusShardDataProvider> {
    // set by the tranversal
    protected int activeRegionExtension = -1;
    protected int maxRegionSize = -1;

    /**
     * our log, which we want to capture anything from this class
     */
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);
    protected final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();

    abstract protected T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker);

    @Override
    public String getTraversalUnits() {
        return "active regions";
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
     * @param activeRegionExtension
     * @param maxRegionSize
     * @return
     */
    protected ActivityProfile incorporateActiveRegions(final ActivityProfile profile,
                                                       final List<ActiveRegion> activeRegions,
                                                       final int activeRegionExtension,
                                                       final int maxRegionSize) {
        if ( profile.isEmpty() )
            throw new IllegalStateException("trying to incorporate an empty active profile " + profile);

        final ActivityProfile bandPassFiltered = profile.bandPassFilter();
        activeRegions.addAll(bandPassFiltered.createActiveRegions( activeRegionExtension, maxRegionSize ));
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

    protected T processActiveRegions(final ActiveRegionWalker<M, T> walker, T sum, final boolean forceRegionsToBeActive) {
        if( walker.activeRegionOutStream != null ) {
            writeActiveRegionsToStream(walker);
            return sum;
        } else {
            return callWalkerMapOnActiveRegions(walker, sum, forceRegionsToBeActive);
        }
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param walker
     */
    private void writeActiveRegionsToStream( final ActiveRegionWalker<M, T> walker ) {
        // Just want to output the active regions to a file, not actually process them
        for( final ActiveRegion activeRegion : workQueue ) {
            if( activeRegion.isActive ) {
                walker.activeRegionOutStream.println( activeRegion.getLocation() );
            }
        }
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
        return (region.getStop() < (getStartOfLiveRegion().getStart() - (includeExtension ? activeRegionExtension : 0)))
                || ! region.onSameContig(getStartOfLiveRegion());
    }

    private T callWalkerMapOnActiveRegions(final ActiveRegionWalker<M, T> walker, T sum, final boolean forceRegionsToBeActive) {
        // Since we've traversed sufficiently past this point (or this contig!) in the workQueue we can unload those regions and process them
        // TODO can implement parallel traversal here
        while( workQueue.peek() != null ) {
            final GenomeLoc extendedLoc = workQueue.peek().getExtendedLoc();
            if ( forceRegionsToBeActive || regionCompletelyWithinDeadZone(extendedLoc, false) ) {
                final ActiveRegion activeRegion = workQueue.remove();
                logger.warn("Processing active region " + activeRegion + " dead zone " + getStartOfLiveRegion());
                sum = processActiveRegion( activeRegion, sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal(final Walker<M, T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M, T>)walker, sum, true);
    }

    protected ActiveRegion getBestRegion(final ActiveRegion activeRegion, final GenomeLoc readLoc) {
        ActiveRegion bestRegion = activeRegion;
        long maxOverlap = activeRegion.getLocation().sizeOfOverlap( readLoc );
        for( final ActiveRegion otherRegionToTest : workQueue ) {
            if( otherRegionToTest.getLocation().sizeOfOverlap(readLoc) >= maxOverlap ) {
                maxOverlap = otherRegionToTest.getLocation().sizeOfOverlap( readLoc );
                bestRegion = otherRegionToTest;
            }
        }
        return bestRegion;
    }
}
