package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
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

    private final Queue<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();
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
        final GenomeLocSortedSet initialIntervals = engine.getIntervals(); // BUGBUG: unfortunate inefficiency that needs to be removed

        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );
        final int activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all

            int minStart = Integer.MAX_VALUE;
            final ArrayList<Double> isActiveList = new ArrayList<Double>();
            GenomeLoc firstIsActiveStart = null;

            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            // We keep processing while the next reference location is within the interval
            GenomeLoc prevLoc = null;
            while( locusView.hasNext() ) {
                final AlignmentContext locus = locusView.next();
                GenomeLoc location = locus.getLocation();
                if(prevLoc != null) {
                    for(int iii = prevLoc.getStart() + 1; iii < location.getStart(); iii++ ) {       
                        final GenomeLoc fakeLoc = engine.getGenomeLocParser().createGenomeLoc(prevLoc.getContig(), iii, iii);
                        if( initialIntervals == null || initialIntervals.overlaps( fakeLoc ) ) {
                            final double isActiveProb = ( walker.presetActiveRegions == null ? 0.0 : ( walker.presetActiveRegions.overlaps(fakeLoc) ? 1.0 : 0.0 ) );
                            isActiveList.add( isActiveProb );
                            if( firstIsActiveStart == null ) {
                                firstIsActiveStart = fakeLoc;
                            }
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
                    final double isActiveProb = ( walker.presetActiveRegions == null ? walker.isActive( tracker, refContext, locus )
                                                                                     : ( walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0 ) );
                    isActiveList.add( isActiveProb );
                    if( firstIsActiveStart == null ) {
                        firstIsActiveStart = location;
                    }
                }

                // Grab all the previously unseen reads from this pileup and add them to the massive read list
                for( final PileupElement p : locus.getBasePileup() ) {
                    final GATKSAMRecord read = p.getRead();
                    if( !myReads.contains(read) ) {
                        myReads.add(read);
                    }
                }

                // If this is the last pileup for this shard calculate the minimum alignment start so that we know 
                // which active regions in the work queue are now safe to process
                if( !locusView.hasNext() ) {
                    for( final PileupElement p : locus.getBasePileup() ) {
                        final GATKSAMRecord read = p.getRead();
                        if( !myReads.contains(read) ) {
                            myReads.add(read);
                        }
                        if( read.getAlignmentStart() < minStart ) { minStart = read.getAlignmentStart(); }
                    }
                }
                prevLoc = location;
                printProgress(dataProvider.getShard(), locus.getLocation());
            }

            // Take the individual isActive calls and integrate them into contiguous active regions and
            // add these blocks of work to the work queue
            final ArrayList<ActiveRegion> activeRegions = integrateActiveList( isActiveList, firstIsActiveStart, activeRegionExtension, walker.presetActiveRegions != null );
            logger.debug("Integrated " + isActiveList.size() + " isActive calls into " + activeRegions.size() + " regions." );
            if( walker.activeRegionOutStream == null ) { 
                workQueue.addAll( activeRegions ); 
            } else { // Just want to output the active regions to a file, not actually process them
                for( final ActiveRegion activeRegion : activeRegions ) {
                    if( activeRegion.isActive ) {
                        walker.activeRegionOutStream.println( activeRegion.getLocation() );
                    }
                }
            }

            // Since we've traversed sufficiently past this point (or this contig!) in the workQueue we can unload those regions and process them
            while( workQueue.peek() != null && (workQueue.peek().getExtendedLoc().getStop() < minStart || !workQueue.peek().getExtendedLoc().getContig().equals(dataProvider.getLocus().getContig())) ) {
                final ActiveRegion activeRegion = workQueue.remove();
                sum = processActiveRegion( activeRegion, myReads, workQueue, sum, walker );
            }
        }

        return sum;
    }

    // Special function called in LinearMicroScheduler to empty out the work queue. Ugly for now but will be cleaned up when we push this functionality more into the engine
    public T endTraversal( final Walker<M,T> walker, T sum) {
        while( workQueue.peek() != null ) {
            final ActiveRegion activeRegion = workQueue.remove();
            sum = processActiveRegion( activeRegion, myReads, workQueue, sum, (ActiveRegionWalker<M,T>) walker );
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

    // band-pass filter the list of isActive probabilities and turn into active regions
    private ArrayList<ActiveRegion> integrateActiveList( final ArrayList<Double> activeList, final GenomeLoc firstIsActiveStart, final int activeRegionExtension, final boolean presetRegions ) {

        final double ACTIVE_PROB_THRESHOLD = 0.2; // BUGBUG: needs to be set-able by the walker author
        final ArrayList<ActiveRegion> returnList = new ArrayList<ActiveRegion>();
        if( activeList.size() == 0 ) {
            return returnList;
        } else if( activeList.size() == 1 ) {
            returnList.add( new ActiveRegion( engine.getGenomeLocParser().createGenomeLoc(firstIsActiveStart.getContig(), firstIsActiveStart.getStart(), firstIsActiveStart.getStart()),
                    activeList.get(0) > ACTIVE_PROB_THRESHOLD, engine.getGenomeLocParser(), activeRegionExtension ) );
            return returnList;
        } else {
            final Double[] activeProbArray = activeList.toArray(new Double[activeList.size()]);
            final double[] filteredProbArray = new double[activeProbArray.length];
            final int FILTER_SIZE = ( presetRegions ? 0 : 50 ); // BUGBUG: needs to be set-able by the walker author
            final int MAX_ACTIVE_REGION = ( presetRegions ? 16001 : 425 ); // BUGBUG: needs to be set-able by the walker author
            for( int iii = 0; iii < activeProbArray.length; iii++ ) {
                double maxVal = 0;
                for( int jjj = Math.max(0, iii-FILTER_SIZE); jjj < Math.min(activeList.size(), iii+FILTER_SIZE+1); jjj++ ) {
                    if( activeProbArray[jjj] > maxVal ) { maxVal = activeProbArray[jjj]; }
                }
                filteredProbArray[iii] = maxVal;
            }
    
            boolean curStatus = filteredProbArray[0] > ACTIVE_PROB_THRESHOLD;
            int curStart = 0;
            for(int iii = 1; iii < filteredProbArray.length; iii++ ) {
                final boolean thisStatus = filteredProbArray[iii] > ACTIVE_PROB_THRESHOLD;
                if( curStatus != thisStatus || (iii-curStart) > MAX_ACTIVE_REGION ) {
                    returnList.add( new ActiveRegion(
                            engine.getGenomeLocParser().createGenomeLoc(firstIsActiveStart.getContig(), firstIsActiveStart.getStart() + curStart, firstIsActiveStart.getStart() + (iii-1)),
                            curStatus, engine.getGenomeLocParser(), activeRegionExtension ) );
                    curStatus = thisStatus;
                    curStart = iii;
                }
            }
            if( curStart != filteredProbArray.length-1 ) {
                returnList.add( new ActiveRegion(
                        engine.getGenomeLocParser().createGenomeLoc(firstIsActiveStart.getContig(), firstIsActiveStart.getStart() + curStart, firstIsActiveStart.getStart() + (filteredProbArray.length-1)),
                        curStatus, engine.getGenomeLocParser(), activeRegionExtension ) );
            }
            return returnList;
        }
    }
}
