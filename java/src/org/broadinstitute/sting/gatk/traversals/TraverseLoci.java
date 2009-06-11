package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.AllLocusView;
import org.broadinstitute.sting.gatk.datasources.providers.CoveredLocusView;
import org.broadinstitute.sting.gatk.datasources.providers.LocusView;
import org.broadinstitute.sting.gatk.datasources.providers.ReferenceOrderedView;
import org.broadinstitute.sting.gatk.datasources.providers.LocusReferenceView;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.ArrayList;
import java.io.File;

/**
 * A simple solution to iterating over all reference positions over a series of genomic locations.
 */
public class TraverseLoci extends TraversalEngine {

    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);

    public TraverseLoci(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super( reads, ref, rods );
    }

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
        logger.debug(String.format("TraverseLoci.traverse Genomic interval is %s", shard.getGenomeLoc()));

        if ( !(walker instanceof LocusWalker) )
            throw new IllegalArgumentException("Walker isn't a loci walker!");

        LocusWalker<M, T> locusWalker = (LocusWalker<M, T>)walker;

        LocusView locusView = getLocusView( walker, dataProvider );
        LocusReferenceView referenceView = new LocusReferenceView( dataProvider );
        ReferenceOrderedView referenceOrderedDataView = new ReferenceOrderedView( dataProvider );

        // We keep processing while the next reference location is within the interval
        while( locusView.hasNext() ) {
            LocusContext locus = locusView.next();

            TraversalStatistics.nRecords++;

            // Iterate forward to get all reference ordered data covering this locus
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation());

            char refBase = referenceView.getReferenceBase(locus.getLocation());

            final boolean keepMeP = locusWalker.filter(tracker, refBase, locus);
            if (keepMeP) {
                M x = locusWalker.map(tracker, refBase, locus);
                sum = locusWalker.reduce(x, sum);
            }

            if (this.maxReads > 0 && TraversalStatistics.nRecords > this.maxReads) {
                logger.warn(String.format("Maximum number of reads encountered, terminating traversal " + TraversalStatistics.nRecords));
                break;
            }

            printProgress("loci", locus.getLocation());
        }

        return sum;
    }

    /**
     * Temporary override of printOnTraversalDone.
     * TODO: Add some sort of TE.getName() function once all TraversalEngines are ported.
     * @param sum Result of the computation.
     * @param <T> Type of the result.
     */
    public <T> void printOnTraversalDone( T sum ) {
        printOnTraversalDone( "loci", sum );
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
        else if( dataSource == DataSource.REFERENCE )
            return new AllLocusView(dataProvider);
        else
            throw new UnsupportedOperationException("Unsupported traversal type: " + dataSource);
    }
}
