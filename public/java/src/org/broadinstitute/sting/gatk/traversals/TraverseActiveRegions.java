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
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
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
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    private final LinkedList<org.broadinstitute.sting.utils.activeregion.ActiveRegion> workQueue = new LinkedList<org.broadinstitute.sting.utils.activeregion.ActiveRegion>();
    private final LinkedHashSet<GATKSAMRecord> myReads = new LinkedHashSet<GATKSAMRecord>();

    @Override
    protected String getTraversalType() {
        return "active regions";
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseActiveRegion.traverse: Shard is %s", dataProvider));

        final LocusView locusView = getLocusView( walker, dataProvider );
        final GenomeLocSortedSet initialIntervals = engine.getIntervals();

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );
        final int activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        final int maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all
            int minStart = Integer.MAX_VALUE;
            ActivityProfile profile = new ActivityProfile(engine.getGenomeLocParser(), walker.hasPresetActiveRegions() );

            ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, dataProvider, locusView);

            // We keep processing while the next reference location is within the interval
            GenomeLoc prevLoc = null;
            while( locusView.hasNext() ) {
                final AlignmentContext locus = locusView.next();
                GenomeLoc location = locus.getLocation();

                if(prevLoc != null) {
                    // fill in the active / inactive labels from the stop of the previous location to the start of this location
                    // TODO refactor to separate function
                    for(int iii = prevLoc.getStop() + 1; iii < location.getStart(); iii++ ) {
                        final GenomeLoc fakeLoc = engine.getGenomeLocParser().createGenomeLoc(prevLoc.getContig(), iii, iii);
                        if( initialIntervals == null || initialIntervals.overlaps( fakeLoc ) ) {
                            final double isActiveProb = ( walker.hasPresetActiveRegions() && walker.presetActiveRegions.overlaps(fakeLoc) ? 1.0 : 0.0 );
                            profile.add(fakeLoc, isActiveProb);
                        }
                    }
                }

                dataProvider.getShard().getReadMetrics().incrementNumIterations();

                // create reference context. Note that if we have a pileup of "extended events", the context will
                // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
                final ReferenceContext refContext = referenceView.getReferenceContext(location);

                // Iterate forward to get all reference ordered data covering this location
                final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

                // Call the walkers isActive function for this locus and add them to the list to be integrated later
                if( initialIntervals == null || initialIntervals.overlaps( location ) ) {
                    final double isActiveProb = walkerActiveProb(walker, tracker, refContext, locus, location);
                    profile.add(location, isActiveProb);
                }

                // Grab all the previously unseen reads from this pileup and add them to the massive read list
                for( final PileupElement p : locus.getBasePileup() ) {
                    final GATKSAMRecord read = p.getRead();
                    if( !myReads.contains(read) ) {
                        myReads.add(read);
                    }

                    // If this is the last pileup for this shard calculate the minimum alignment start so that we know
                    // which active regions in the work queue are now safe to process
                    minStart = Math.min(minStart, read.getAlignmentStart());
                }

                prevLoc = location;

                printProgress(dataProvider.getShard(), locus.getLocation());
            }

            // Take the individual isActive calls and integrate them into contiguous active regions and
            // add these blocks of work to the work queue
            // band-pass filter the list of isActive probabilities and turn into active regions
            final ActivityProfile bandPassFiltered = profile.bandPassFilter();
            final List<org.broadinstitute.sting.utils.activeregion.ActiveRegion> activeRegions = bandPassFiltered.createActiveRegions( activeRegionExtension, maxRegionSize );

            // add active regions to queue of regions to process
            // first check if can merge active regions over shard boundaries
            if( !activeRegions.isEmpty() ) {
                if( !workQueue.isEmpty() ) {
                    final org.broadinstitute.sting.utils.activeregion.ActiveRegion last = workQueue.getLast();
                    final org.broadinstitute.sting.utils.activeregion.ActiveRegion first = activeRegions.get(0);
                    if( last.isActive == first.isActive && last.getLocation().contiguousP(first.getLocation()) && last.getLocation().size() + first.getLocation().size() <= maxRegionSize ) {
                        workQueue.removeLast();
                        activeRegions.remove(first);
                        workQueue.add( new org.broadinstitute.sting.utils.activeregion.ActiveRegion(last.getLocation().union(first.getLocation()), first.isActive, this.engine.getGenomeLocParser(), activeRegionExtension) );
                    }
                }
                workQueue.addAll( activeRegions );
            }

            logger.debug("Integrated " + profile.size() + " isActive calls into " + activeRegions.size() + " regions." );

            // now go and process all of the active regions
            sum = processActiveRegions(walker, sum, minStart, dataProvider.getLocus().getContig());
        }

        return sum;
    }


    // --------------------------------------------------------------------------------
    //
    // simple utility functions
    //
    // --------------------------------------------------------------------------------

    private final double walkerActiveProb(final ActiveRegionWalker<M,T> walker,
                                          final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                          final AlignmentContext locus, final GenomeLoc location) {
        if ( walker.hasPresetActiveRegions() ) {
            return walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0;
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
        for( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion : workQueue ) {
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
                final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion = workQueue.remove();
                sum = processActiveRegion( activeRegion, myReads, workQueue, sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    private T processActiveRegion( final org.broadinstitute.sting.utils.activeregion.ActiveRegion activeRegion, final LinkedHashSet<GATKSAMRecord> reads, final Queue<org.broadinstitute.sting.utils.activeregion.ActiveRegion> workQueue, final T sum, final ActiveRegionWalker<M,T> walker ) {
        final ArrayList<GATKSAMRecord> placedReads = new ArrayList<GATKSAMRecord>();
        for( final GATKSAMRecord read : reads ) {
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );
            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                // The region which the highest amount of overlap is chosen as the primary region for the read (tie breaking is done as right most region)
                long maxOverlap = activeRegion.getLocation().sizeOfOverlap( readLoc );
                org.broadinstitute.sting.utils.activeregion.ActiveRegion bestRegion = activeRegion;
                for( final org.broadinstitute.sting.utils.activeregion.ActiveRegion otherRegionToTest : workQueue ) {
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
                    for( final org.broadinstitute.sting.utils.activeregion.ActiveRegion otherRegionToTest : workQueue ) {
                        if( !bestRegion.equals(otherRegionToTest) && otherRegionToTest.getExtendedLoc().overlapsP( readLoc ) ) {
                            otherRegionToTest.add( read );
                        }
                    }
                }
                placedReads.add( read );
            } else if( activeRegion.getExtendedLoc().overlapsP( readLoc ) && walker.wantsNonPrimaryReads() ) {
                activeRegion.add( read );
            }
        }
        reads.removeAll( placedReads ); // remove all the reads which have been placed into their active region

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLoc());
        final M x = walker.map( activeRegion, null );
        return walker.reduce( x, sum );
    }

    // --------------------------------------------------------------------------------
    //
    // engine interaction code
    //
    // --------------------------------------------------------------------------------

    /**
     * Gets the best view of loci for this walker given the available data.
     * @param walker walker to interrogate.
     * @param dataProvider Data which which to drive the locus view.
     * @return A view of the locus data, where one iteration of the locus view maps to one iteration of the traversal.
     */
    private LocusView getLocusView( final Walker<M,T> walker, final LocusShardDataProvider dataProvider ) {
        final DataSource dataSource = WalkerManager.getWalkerDataSource(walker);
        if( dataSource == DataSource.READS )
            return new CoveredLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE ) //|| ! GenomeAnalysisEngine.instance.getArguments().enableRodWalkers )
            return new AllLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE_ORDERED_DATA )
            return new RodLocusView(dataProvider);
        else
            throw new UnsupportedOperationException("Unsupported traversal type: " + dataSource);
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal( final Walker<M,T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M,T>)walker, sum, Integer.MAX_VALUE, null);
    }
}
