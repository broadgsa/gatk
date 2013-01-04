package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.executive.WindowMaker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

public class ExperimentalActiveRegionShardTraverseActiveRegions <M,T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,ActiveRegionShardDataProvider> {
    /**
     * our log, which we want to capture anything from this class
     */
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);

    private final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();
    private final LinkedList<GATKSAMRecord> myReads = new LinkedList<GATKSAMRecord>();

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final ActiveRegionShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("ExperimentalActiveRegionShardTraverseActiveRegions.traverse: Shard is %s", dataProvider));

        ReadShardDataProvider readDataProvider = dataProvider.getReadShardDataProvider();

        final int activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        final int maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

        final ReadView readView = new ReadView(readDataProvider);

        final List<ActiveRegion> activeRegions = new LinkedList<ActiveRegion>();
        ActivityProfile profile = new ActivityProfile(engine.getGenomeLocParser(), walker.hasPresetActiveRegions());

        Shard readShard = readDataProvider.getShard();
        SAMFileHeader header = readShard.getReadProperties().getHeader();
        WindowMaker windowMaker = new WindowMaker(readShard, engine.getGenomeLocParser(),
                readView.iterator(), readShard.getGenomeLocs(), SampleUtils.getSAMFileSamples(header));

        for(WindowMaker.WindowMakerIterator iterator: windowMaker) {
            LocusShardDataProvider locusDataProvider = dataProvider.getLocusShardDataProvider(iterator);
            final LocusView locusView = new AllLocusView(locusDataProvider);
            final LocusReferenceView referenceView = new LocusReferenceView( walker, locusDataProvider );
            ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, locusDataProvider, locusView);

            // We keep processing while the next reference location is within the interval
            GenomeLoc prevLoc = null;
            while( locusView.hasNext() ) {
                final AlignmentContext locus = locusView.next();
                final GenomeLoc location = locus.getLocation();

                if ( prevLoc != null && location.getStart() != prevLoc.getStop() + 1 ) {
                    // we've move across some interval boundary, restart profile
                    profile = incorporateActiveRegions(profile, activeRegions, activeRegionExtension, maxRegionSize);
                }

                readDataProvider.getShard().getReadMetrics().incrementNumIterations();

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

            locusDataProvider.close();
        }

        windowMaker.close();

        updateCumulativeMetrics(readDataProvider.getShard());

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
                    workQueue.addLast(new ActiveRegion(last.getLocation().union(first.getLocation()), first.isActive, this.engine.getGenomeLocParser(), activeRegionExtension));
                }
            }
            workQueue.addAll( activeRegions );
        }

        logger.debug("Integrated " + profile.size() + " isActive calls into " + activeRegions.size() + " regions." );

        // now process the active regions, where possible
        boolean emptyQueue = false;
        sum = processActiveRegions(walker, sum, emptyQueue);

        return sum;
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

    private T processActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, boolean emptyQueue ) {
        if( walker.activeRegionOutStream != null ) {
            writeActiveRegionsToStream(walker);
            return sum;
        } else {
            return callWalkerMapOnActiveRegions(walker, sum, emptyQueue);
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

    private T callWalkerMapOnActiveRegions( final ActiveRegionWalker<M,T> walker, T sum, boolean emptyQueue ) {
        final int lastRegionStart = workQueue.getLast().getLocation().getStart();
        final String lastRegionContig = workQueue.getLast().getLocation().getContig();

        // If we've traversed sufficiently past the beginning of the workQueue we can unload those regions and process them
        // TODO can implement parallel traversal here
        while( workQueue.peekFirst() != null ) {
            ActiveRegion firstRegion = workQueue.getFirst();
            final String firstRegionContig = firstRegion.getLocation().getContig();
            if (emptyQueue || firstRegionContig != lastRegionContig) {
                sum = processFirstActiveRegion(sum, walker);
            }
            else {
                final int firstRegionMaxReadStop = walker.wantsExtendedReads() ? firstRegion.getMaxReadStop() : firstRegion.getExtendedMaxReadStop();
                if (lastRegionStart > firstRegionMaxReadStop) {
                    sum = processFirstActiveRegion( sum, walker );
                }
                else {
                    break;
                }
            }
        }

        return sum;
    }

    /**
     * Process the first active region and all remaining reads which overlap
     *
     * Remove the first active region from the queue
     * (NB: some reads associated with this active region may have already been processed)
     *
     * Remove all of these reads from the queue
     * (NB: some may be associated with other active regions)
     *
     * @param sum
     * @param walker
     * @return
     */
    private T processFirstActiveRegion( final T sum, final ActiveRegionWalker<M,T> walker ) {
        final ActiveRegion firstRegion = workQueue.removeFirst();

        GATKSAMRecord firstRead = myReads.peekFirst();     // don't remove because it may not be placed here
        GenomeLoc firstReadLoc = this.engine.getGenomeLocParser().createGenomeLoc( firstRead );

        while ( firstRegion.getLocation().overlapsP( firstReadLoc ) ||
                (walker.wantsExtendedReads() && firstRegion.getExtendedLoc().overlapsP( firstReadLoc ))) {
            if( firstRegion.getLocation().overlapsP( firstReadLoc ) ) {
                // The region which the highest amount of overlap is chosen as the primary region for the read (tie breaking is done as right most region)
                long maxOverlap = firstRegion.getLocation().sizeOfOverlap( firstReadLoc );
                ActiveRegion bestRegion = firstRegion;
                for( final ActiveRegion otherRegionToTest : workQueue ) {
                    if( otherRegionToTest.getLocation().sizeOfOverlap(firstReadLoc) >= maxOverlap ) {
                        maxOverlap = otherRegionToTest.getLocation().sizeOfOverlap( firstReadLoc );
                        bestRegion = otherRegionToTest;
                    }
                }
                bestRegion.add( firstRead );

                // The read is also added to all other regions in which it overlaps but marked as non-primary
                if( walker.wantsNonPrimaryReads() ) {
                    if( !bestRegion.equals(firstRegion) ) {
                        firstRegion.add(firstRead);
                    }
                    for( final ActiveRegion otherRegionToTest : workQueue ) {
                        if( !bestRegion.equals(otherRegionToTest) ) {
                            // check for non-primary vs. extended
                            if ( otherRegionToTest.getLocation().overlapsP( firstReadLoc ) ) {
                                otherRegionToTest.add( firstRead );
                            } else if ( walker.wantsExtendedReads() && otherRegionToTest.getExtendedLoc().overlapsP( firstReadLoc ) ) {
                                otherRegionToTest.add( firstRead );
                            }
                        }
                    }
                }

                // check for non-primary vs. extended
            } else if( firstRegion.getLocation().overlapsP( firstReadLoc ) ) {
                if ( walker.wantsNonPrimaryReads() ) {
                    firstRegion.add( firstRead );
                }
            } else if( walker.wantsExtendedReads() && firstRegion.getExtendedLoc().overlapsP( firstReadLoc )) {
                firstRegion.add( firstRead );
            }

            myReads.removeFirst();
            firstRead = myReads.peekFirst();
            firstReadLoc = this.engine.getGenomeLocParser().createGenomeLoc( firstRead );
        }

        logger.debug(">> Map call with " + firstRegion.getReads().size() + " " + (firstRegion.isActive ? "active" : "inactive") + " reads @ " + firstRegion.getLocation() + " with full extent: " + firstRegion.getReferenceLoc());
        final M x = walker.map( firstRegion, null );
        return walker.reduce(x, sum);
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal( final Walker<M,T> walker, T sum) {
        boolean emptyQueue = true;
        return processActiveRegions((ActiveRegionWalker<M,T>)walker, sum, emptyQueue);
    }
}
