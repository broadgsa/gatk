package org.broadinstitute.sting.gatk.traversals;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Queue;

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
    private final LinkedHashSet<SAMRecord> myReads = new LinkedHashSet<SAMRecord>();

    @Override
    protected String getTraversalType() {
        return "active regions";
    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseActiveRegion.traverse: Shard is %s", dataProvider));

        LocusView locusView = getLocusView( walker, dataProvider );

        int minStart = Integer.MAX_VALUE;
        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all

            final ArrayList<ActiveRegion> isActiveList = new ArrayList<ActiveRegion>();

            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            // We keep processing while the next reference location is within the interval
            while( locusView.hasNext() ) {
                final AlignmentContext locus = locusView.next();
                GenomeLoc location = locus.getLocation();

                dataProvider.getShard().getReadMetrics().incrementNumIterations();

                if ( locus.hasExtendedEventPileup() ) {
                    // if the alignment context we received holds an "extended" pileup (i.e. pileup of insertions/deletions
                    // associated with the current site), we need to update the location. The updated location still starts
                    // at the current genomic position, but it has to span the length of the longest deletion (if any).
                    location = engine.getGenomeLocParser().setStop(location,location.getStop()+locus.getExtendedEventPileup().getMaxDeletionLength());

                    // it is possible that the new expanded location spans the current shard boundary; the next method ensures
                    // that when it is the case, the reference sequence held by the ReferenceView will be reloaded so that
                    // the view has all the bases we are gonna need. If the location fits within the current view bounds,
                    // the next call will not do anything to the view:
                    referenceView.expandBoundsToAccomodateLoc(location);
                }

                // create reference context. Note that if we have a pileup of "extended events", the context will
                // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
                final ReferenceContext refContext = referenceView.getReferenceContext(location);

                // Iterate forward to get all reference ordered data covering this location
                final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

                // Call the walkers isActive function for this locus and add them to the list to be integrated later
                final boolean isActive = walker.isActive( tracker, refContext, locus );
                isActiveList.add( new ActiveRegion(location, isActive, engine.getGenomeLocParser()) );

                // Grab all the previously unseen reads from this pileup and add them to the massive read list
                for( final PileupElement p : locus.getBasePileup() ) {
                    final SAMRecord read = p.getRead();
                    if( !myReads.contains(read) ) {
                        myReads.add(read);
                    }
                }

                // If this is the last pileup for this shard then need to calculate the minimum alignment start so that
                // we know which active regions in the work queue are now safe to process
                if( !locusView.hasNext() ) {
                    for( final PileupElement p : locus.getBasePileup() ) {
                        final SAMRecord read = p.getRead();
                        if( read.getAlignmentStart() < minStart ) { minStart = read.getAlignmentStart(); }
                    }
                }
                printProgress(dataProvider.getShard(),locus.getLocation());
            }

            // Take the individual isActive calls and integrate them into contiguous active regions and
            // add these blocks of work to the work queue
            final ArrayList<ActiveRegion> activeRegions = integrateActiveList( isActiveList );
            logger.debug("Integrated " + isActiveList.size() + " isActive calls into " + activeRegions.size() + " regions." );
            workQueue.addAll( activeRegions );
        }

        while( workQueue.peek().getLocation().getStop() < minStart ) {
            final ActiveRegion activeRegion = workQueue.remove();
            sum = processActiveRegion( activeRegion, myReads, workQueue, sum, walker );
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

    private T processActiveRegion( final ActiveRegion activeRegion, final LinkedHashSet<SAMRecord> reads, final Queue<ActiveRegion> workQueue, final T sum, final ActiveRegionWalker<M,T> walker ) {
        final ArrayList<SAMRecord> placedReads = new ArrayList<SAMRecord>();
        for( final SAMRecord read : reads ) {
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
                bestRegion.add( (GATKSAMRecord) read, true );

                // The read is also added to all other regions in which it overlaps but marked as non-primary
                if( !bestRegion.equals(activeRegion) ) {
                    activeRegion.add( (GATKSAMRecord) read, false );
                }
                for( final ActiveRegion otherRegionToTest : workQueue ) {
                    if( !bestRegion.equals(otherRegionToTest) && otherRegionToTest.getLocation().overlapsP( readLoc ) ) {
                        activeRegion.add( (GATKSAMRecord) read, false );
                    }
                }
                placedReads.add( read );
            }
        }
        reads.removeAll( placedReads ); // remove all the reads which have been placed into their active region

        logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReferenceLocation());
        final M x = walker.map( activeRegion, null ); // BUGBUG: tracker needs to be filled in and passed to the walker
        return walker.reduce( x, sum );
    }

    /**
     * Gets the best view of loci for this walker given the available data.
     * @param walker walker to interrogate.
     * @param dataProvider Data which which to drive the locus view.
     * @return A view of the locus data, where one iteration of the locus view maps to one iteration of the traversal.
     */
    private LocusView getLocusView( Walker<M,T> walker, LocusShardDataProvider dataProvider ) {
        DataSource dataSource = WalkerManager.getWalkerDataSource(walker);
        if( dataSource == DataSource.READS )
            return new CoveredLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE ) //|| ! GenomeAnalysisEngine.instance.getArguments().enableRodWalkers )
            return new AllLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE_ORDERED_DATA )
            return new RodLocusView(dataProvider);
        else
            throw new UnsupportedOperationException("Unsupported traversal type: " + dataSource);
    }

    // integrate active regions into contiguous chunks based on active status
    private ArrayList<ActiveRegion> integrateActiveList( final ArrayList<ActiveRegion> activeList ) {
        final ArrayList<ActiveRegion> returnList = new ArrayList<ActiveRegion>();
        ActiveRegion prevLocus = activeList.remove(0);
        ActiveRegion startLocus = prevLocus;
        for( final ActiveRegion thisLocus : activeList ) {
            if( prevLocus.isActive != thisLocus.isActive ) {
                returnList.add( new ActiveRegion( engine.getGenomeLocParser().createGenomeLoc(startLocus.getLocation().getContig(), startLocus.getLocation().getStart(), prevLocus.getLocation().getStart()),
                                                  prevLocus.isActive, engine.getGenomeLocParser() ) );
                startLocus = thisLocus;
            }
            prevLocus = thisLocus;
        }
        // output the last region if necessary
        if( startLocus != prevLocus ) {
            returnList.add( new ActiveRegion( engine.getGenomeLocParser().createGenomeLoc(startLocus.getLocation().getContig(), startLocus.getLocation().getStart(), prevLocus.getLocation().getStart()),
                                              prevLocus.isActive, engine.getGenomeLocParser() ) );
        }
        return returnList;
    }
}
