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
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 12/9/11
 */

public class TraverseActiveRegions <M,T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,LocusShardDataProvider> {
    /**
     * our log, which we want to capture anything from this class
     */
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);

    private final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();
    private final LinkedHashSet<GATKSAMRecord> myReads = new LinkedHashSet<GATKSAMRecord>();

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseActiveRegion.traverse: Shard is %s", dataProvider));

        final LocusView locusView = new AllLocusView(dataProvider);

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );
        final int activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        final int maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

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
                profile = incorporateActiveRegions(profile, activeRegions, activeRegionExtension, maxRegionSize);
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
            incorporateActiveRegions(profile, activeRegions, activeRegionExtension, maxRegionSize);

        // add active regions to queue of regions to process
        // first check if can merge active regions over shard boundaries
        if( !activeRegions.isEmpty() ) {
            if( !workQueue.isEmpty() ) {
                final ActiveRegion last = workQueue.getLast();
                final ActiveRegion first = activeRegions.get(0);
                if( last.isActive == first.isActive && last.getLocation().contiguousP(first.getLocation()) && last.getLocation().size() + first.getLocation().size() <= maxRegionSize ) {
                    workQueue.removeLast();
                    activeRegions.remove(first);
                    workQueue.add( new ActiveRegion(last.getLocation().union(first.getLocation()), first.isActive, this.engine.getGenomeLocParser(), activeRegionExtension) );
                }
            }
            workQueue.addAll( activeRegions );
        }

        logger.debug("Integrated " + profile.size() + " isActive calls into " + activeRegions.size() + " regions." );

        // now go and process all of the active regions
        sum = processActiveRegions(walker, sum, minStart, dataProvider.getLocus().getContig());

        return sum;
    }

    /**
     * Is the loc outside of the intervals being requested for processing by the GATK?
     * @param loc
     * @return
     */
    private boolean outsideEngineIntervals(final GenomeLoc loc) {
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
    private ActivityProfile incorporateActiveRegions(final ActivityProfile profile,
                                                     final List<ActiveRegion> activeRegions,
                                                     final int activeRegionExtension,
                                                     final int maxRegionSize) {
        if ( profile.isEmpty() )
            throw new IllegalStateException("trying to incorporate an empty active profile " + profile);

        final ActivityProfile bandPassFiltered = profile.bandPassFilter();
        activeRegions.addAll(bandPassFiltered.createActiveRegions( activeRegionExtension, maxRegionSize ));
        return new ActivityProfile( engine.getGenomeLocParser(), profile.hasPresetRegions() );
    }


    // --------------------------------------------------------------------------------
    //
    // simple utility functions
    //
    // --------------------------------------------------------------------------------

    private final ActivityProfileResult walkerActiveProb(final ActiveRegionWalker<M,T> walker,
                                          final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                          final AlignmentContext locus, final GenomeLoc location) {
        if ( walker.hasPresetActiveRegions() ) {
            return new ActivityProfileResult(location, walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0);
        } else {
            return walker.isActive( tracker, refContext, locus );
        }
    }

    private ReferenceOrderedView getReferenceOrderedView( final ActiveRegionWalker<M,T> walker,
                                                          final LocusShardDataProvider dataProvider,
                                                          final LocusView locusView) {
        if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
            return new ManagingReferenceOrderedView( dataProvider );
        else
            return (RodLocusView)locusView;
    }

    // --------------------------------------------------------------------------------
    //
    // code to handle processing active regions
    //
    // --------------------------------------------------------------------------------

    private T processActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, final int minStart, final String currentContig ) {
        if( walker.activeRegionOutStream != null ) {
            writeActiveRegionsToStream(walker);
            return sum;
        } else {
            return callWalkerMapOnActiveRegions(walker, sum, minStart, currentContig);
        }
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param walker
     */
    private void writeActiveRegionsToStream( final ActiveRegionWalker<M,T> walker ) {
        // Just want to output the active regions to a file, not actually process them
        for( final ActiveRegion activeRegion : workQueue ) {
            if( activeRegion.isActive ) {
                walker.activeRegionOutStream.println( activeRegion.getLocation() );
            }
        }
    }

    private T callWalkerMapOnActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, final int minStart, final String currentContig ) {
        // Since we've traversed sufficiently past this point (or this contig!) in the workQueue we can unload those regions and process them
        // TODO can implement parallel traversal here
        while( workQueue.peek() != null ) {
            final GenomeLoc extendedLoc = workQueue.peek().getExtendedLoc();
            if ( extendedLoc.getStop() < minStart || (currentContig != null && !workQueue.peek().getExtendedLoc().getContig().equals(currentContig))) {
                final ActiveRegion activeRegion = workQueue.remove();
                sum = processActiveRegion( activeRegion, myReads, workQueue, sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    private T processActiveRegion( final ActiveRegion activeRegion, final LinkedHashSet<GATKSAMRecord> reads, final Queue<ActiveRegion> workQueue, final T sum, final ActiveRegionWalker<M,T> walker ) {
        final ArrayList<GATKSAMRecord> placedReads = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord read : reads ) {
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );
            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                // The region which the highest amount of overlap is chosen as the primary region for the read (tie breaking is done as right most region)
                long maxOverlap = activeRegion.getLocation().sizeOfOverlap( readLoc );
                ActiveRegion bestRegion = activeRegion;
                for( final ActiveRegion otherRegionToTest : workQueue ) {
                    if( otherRegionToTest.getLocation().sizeOfOverlap(readLoc) >= maxOverlap ) {
                        maxOverlap = otherRegionToTest.getLocation().sizeOfOverlap( readLoc );
                        bestRegion = otherRegionToTest;
                    }
                }
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
            } else if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                if ( walker.wantsNonPrimaryReads() ) {
                    activeRegion.add( read );
                }
            } else if( walker.wantsExtendedReads() && activeRegion.getExtendedLoc().overlapsP( readLoc )) {
                activeRegion.add( read );
            }
        }
        reads.removeAll( placedReads ); // remove all the reads which have been placed into their active region
        // WARNING: This hashset relies on reads being exactly equal when they are placed in the list as when they are removed. So the ActiveRegionWalker can't modify the reads in any way.

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLoc());
        final M x = walker.map( activeRegion, null );
        return walker.reduce( x, sum );
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal( final Walker<M,T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M,T>)walker, sum, Integer.MAX_VALUE, null);
    }
}
