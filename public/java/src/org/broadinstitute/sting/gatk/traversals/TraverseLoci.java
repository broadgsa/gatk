package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLoci<M,T> extends TraversalEngine<M,T,LocusWalker<M,T>,LocusShardDataProvider> {
    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    @Override
    protected String getTraversalType() {
        return "sites";
    }

    @Override
    public T traverse( LocusWalker<M,T> walker,
                       LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseLoci.traverse: Shard is %s", dataProvider));

        LocusView locusView = getLocusView( walker, dataProvider );
        boolean done = false;

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all

            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

            // We keep processing while the next reference location is within the interval
            while( locusView.hasNext() && ! done ) {
                AlignmentContext locus = locusView.next();
                GenomeLoc location = locus.getLocation();

                dataProvider.getShard().getReadMetrics().incrementNumIterations();

                // create reference context. Note that if we have a pileup of "extended events", the context will
                // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
                ReferenceContext refContext = referenceView.getReferenceContext(location);

                // Iterate forward to get all reference ordered data covering this location
                final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

                final boolean keepMeP = walker.filter(tracker, refContext, locus);
                if (keepMeP) {
                    M x = walker.map(tracker, refContext, locus);
                    sum = walker.reduce(x, sum);
                    done = walker.isDone();
                }

                printProgress(dataProvider.getShard(),locus.getLocation());
            }
        }

        // We have a final map call to execute here to clean up the skipped based from the
        // last position in the ROD to that in the interval
        if ( WalkerManager.getWalkerDataSource(walker) == DataSource.REFERENCE_ORDERED_DATA && ! walker.isDone() ) {
            // only do this if the walker isn't done!
            RodLocusView rodLocusView = (RodLocusView)locusView;
            long nSkipped = rodLocusView.getLastSkippedBases();
            if ( nSkipped > 0 ) {
                GenomeLoc site = rodLocusView.getLocOneBeyondShard();
                AlignmentContext ac = new AlignmentContext(site, new ReadBackedPileupImpl(site), nSkipped);
                M x = walker.map(null, null, ac);
                sum = walker.reduce(x, sum);
            }
        }

        return sum;
    }

    /**
     * Gets the best view of loci for this walker given the available data.  The view will function as a 'trigger track'
     * of sorts, providing a consistent interface so that TraverseLoci doesn't need to be reimplemented for any new datatype
     * that comes along.
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
}
