package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLoci extends TraversalEngine {
    final private static String LOCI_STRING = "sites";
    final private static boolean ENABLE_ROD_TRAVERSAL = false;


    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    public <M,T> T traverse(Walker<M,T> walker, ArrayList<GenomeLoc> locations) {
        if ( locations.isEmpty() )
            Utils.scareUser("Requested all locations be processed without providing locations to be processed!");

        throw new UnsupportedOperationException("This traversal type not supported by TraverseLoci");
    }

    @Override
    public <M,T> T traverse( Walker<M,T> walker,
                             Shard shard,
                             ShardDataProvider dataProvider,
                             T sum ) {
        logger.debug(String.format("TraverseLoci.traverse: Shard is %s", shard));

        if ( !(walker instanceof LocusWalker) )
            throw new IllegalArgumentException("Walker isn't a loci walker!");

        LocusWalker<M, T> locusWalker = (LocusWalker<M, T>)walker;

        LocusView locusView = getLocusView( walker, dataProvider );

        //if ( WalkerManager.getWalkerDataSource(walker) == DataSource.REFERENCE_ORDERED_DATA )
        //    throw new RuntimeException("Engine currently doesn't support RodWalkers");

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all

            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( ! ENABLE_ROD_TRAVERSAL || WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

            // We keep processing while the next reference location is within the interval
            while( locusView.hasNext() ) {
                AlignmentContext locus = locusView.next();

                TraversalStatistics.nRecords++;

                // Iterate forward to get all reference ordered data covering this locus
                final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation());

                ReferenceContext refContext = referenceView.getReferenceContext(locus.getLocation());

                final boolean keepMeP = locusWalker.filter(tracker, refContext, locus);
                if (keepMeP) {
                    M x = locusWalker.map(tracker, refContext, locus);
                    sum = locusWalker.reduce(x, sum);
                }

                if (this.maximumIterations > 0 && TraversalStatistics.nRecords > this.maximumIterations) {
                    logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                    break;
                }

                printProgress(LOCI_STRING, locus.getLocation());
            }
        }

            // We have a final map call to execute here to clean up the skipped based from the
            // last position in the ROD to that in the interval
        if ( ENABLE_ROD_TRAVERSAL && WalkerManager.getWalkerDataSource(walker) == DataSource.REFERENCE_ORDERED_DATA ) {
            RodLocusView rodLocusView = (RodLocusView)locusView;
            long nSkipped = rodLocusView.getLastSkippedBases();
            if ( nSkipped > 0 ) {
                // no sense in making the call if you don't have anything interesting to say
                AlignmentContext ac = new AlignmentContext(rodLocusView.getLocOneBeyondShard(), null, null, nSkipped);
                M x = locusWalker.map(null, null, ac);
                sum = locusWalker.reduce(x, sum);
            }
        }

        return sum;
    }

    /**
     * Temporary override of printOnTraversalDone.
     * 
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone(LOCI_STRING, sum );
    }

    /**
     * Gets the best view of loci for this walker given the available data.
     * @param walker walker to interrogate.
     * @param dataProvider Data which which to drive the locus view.
     */
    private LocusView getLocusView( Walker walker, ShardDataProvider dataProvider ) {
        DataSource dataSource = WalkerManager.getWalkerDataSource(walker);
        if( dataSource == DataSource.READS )
            return new CoveredLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE || ! ENABLE_ROD_TRAVERSAL )
            return new AllLocusView(dataProvider);
        else if( dataSource == DataSource.REFERENCE_ORDERED_DATA )
            return new RodLocusView(dataProvider);
        else
            throw new UnsupportedOperationException("Unsupported traversal type: " + dataSource);
    }
}
