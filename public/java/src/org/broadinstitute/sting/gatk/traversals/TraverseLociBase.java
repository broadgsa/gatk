package org.broadinstitute.sting.gatk.traversals;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public abstract class TraverseLociBase<M,T> extends TraversalEngine<M,T,LocusWalker<M,T>,LocusShardDataProvider> {
    /**
     * our log, which we want to capture anything from this class
     */
    protected static final Logger logger = Logger.getLogger(TraversalEngine.class);

    @Override
    protected final String getTraversalType() {
        return "sites";
    }

    protected static class TraverseResults<T> {
        final int numIterations;
        final T reduceResult;

        public TraverseResults(int numIterations, T reduceResult) {
            this.numIterations = numIterations;
            this.reduceResult = reduceResult;
        }
    }

    protected abstract TraverseResults<T> traverse( final LocusWalker<M,T> walker,
                                                    final LocusView locusView,
                                                    final LocusReferenceView referenceView,
                                                    final ReferenceOrderedView referenceOrderedDataView,
                                                    final T sum);

    @Override
    public T traverse( LocusWalker<M,T> walker,
                       LocusShardDataProvider dataProvider,
                       T sum) {
        logger.debug(String.format("TraverseLociBase.traverse: Shard is %s", dataProvider));

        final LocusView locusView = getLocusView( walker, dataProvider );

        if ( locusView.hasNext() ) { // trivial optimization to avoid unnecessary processing when there's nothing here at all
            //ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );
            ReferenceOrderedView referenceOrderedDataView = null;
            if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
                referenceOrderedDataView = new ManagingReferenceOrderedView( dataProvider );
            else
                referenceOrderedDataView = (RodLocusView)locusView;

            final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );

            final TraverseResults<T> result = traverse( walker, locusView, referenceView, referenceOrderedDataView, sum );
            sum = result.reduceResult;
            dataProvider.getShard().getReadMetrics().incrementNumIterations(result.numIterations);
        }

        // We have a final map call to execute here to clean up the skipped based from the
        // last position in the ROD to that in the interval
        if ( WalkerManager.getWalkerDataSource(walker) == DataSource.REFERENCE_ORDERED_DATA && ! walker.isDone() ) {
            // only do this if the walker isn't done!
            final RodLocusView rodLocusView = (RodLocusView)locusView;
            final long nSkipped = rodLocusView.getLastSkippedBases();
            if ( nSkipped > 0 ) {
                final GenomeLoc site = rodLocusView.getLocOneBeyondShard();
                final AlignmentContext ac = new AlignmentContext(site, new ReadBackedPileupImpl(site), nSkipped);
                final M x = walker.map(null, null, ac);
                sum = walker.reduce(x, sum);
            }
        }

        return sum;
    }

    /**
     * Gets the best view of loci for this walker given the available data.  The view will function as a 'trigger track'
     * of sorts, providing a consistent interface so that TraverseLociBase doesn't need to be reimplemented for any new datatype
     * that comes along.
     * @param walker walker to interrogate.
     * @param dataProvider Data which which to drive the locus view.
     * @return A view of the locus data, where one iteration of the locus view maps to one iteration of the traversal.
     */
    private LocusView getLocusView( Walker<M,T> walker, LocusShardDataProvider dataProvider ) {
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
}
