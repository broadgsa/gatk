package org.broadinstitute.sting.gatk.traversals;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.InvalidPositionException;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.ArrayList;
import java.io.File;

/**
 * A simple, short-term solution to iterating over all reference positions over a series of
 * genomic locations. Simply overloads the superclass traverse function to go over the entire
 * interval's reference positions.
 * mhanna - Added better data source integration.
 * TODO: Gain confidence in this implementation and remove the original.
 */
public class TraverseLociByReference extends TraversalEngine {

    /**
     * our log, which we want to capture anything from this class
     */
    protected static Logger logger = Logger.getLogger(TraversalEngine.class);


    public TraverseLociByReference(List<File> reads, File ref, List<ReferenceOrderedData<? extends ReferenceOrderedDatum>> rods) {
        super( reads, ref, rods );
    }

    public <M,T> T traverse(Walker<M,T> walker, ArrayList<GenomeLoc> locations) {
        if ( locations.isEmpty() )
            Utils.scareUser("Requested all locations be processed without providing locations to be processed!");

        throw new UnsupportedOperationException("This traversal type not supported by TraverseLociByReference");
    }

    public <M,T> T traverse( Walker<M,T> walker,
                             Shard shard,
                             ReferenceProvider referenceProvider,
                             LocusContextProvider locusProvider,
                             T sum ) {
        logger.debug(String.format("TraverseLociByReference.traverse Genomic interval is %s", shard.getGenomeLoc()));

        if ( !(walker instanceof LocusWalker) )
            throw new IllegalArgumentException("Walker isn't a loci walker!");

        LocusWalker<M, T> locusWalker = (LocusWalker<M, T>)walker;
        GenomeLoc loc = shard.getGenomeLoc();

        // We keep processing while the next reference location is within the interval
        for( long pos = loc.getStart(); pos <= loc.getStop(); pos++ ) {
            GenomeLoc site = new GenomeLoc( loc.getContig(), pos );

            TraversalStatistics.nRecords++;

            // Iterate forward to get all reference ordered data covering this locus
            final RefMetaDataTracker tracker = getReferenceOrderedDataAtLocus( site );

            LocusContext locus = locusProvider.getLocusContext( site );

            char refBase = referenceProvider.getReferenceBase( site );

            if ( DOWNSAMPLE_BY_COVERAGE )
                locus.downsampleToCoverage(downsamplingCoverage);

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
}
