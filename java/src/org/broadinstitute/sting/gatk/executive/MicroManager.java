package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

import java.util.List;
import java.util.Arrays;
import java.util.Iterator;
import java.io.File;
import java.io.IOException;

/**
 * A micro-scheduling manager for N-way threaded execution of a traversal
 *
 */
public class MicroManager {
    private static long SHARD_SIZE = 5L;

    private File reads;
    private FastaSequenceFile2 ref;

    private TraverseLociByReference traversalEngine = null;

    protected static Logger logger = Logger.getLogger(MicroManager.class);

    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    public MicroManager( File reads,                          // the reads file
                         File refFile,                        // the reference file driving the traversal
                         int nThreadsToUse ) {                // maximum number of threads to use to do the work

        this.reads = reads;
        ref = new FastaSequenceFile2(refFile);
        GenomeLoc.setupRefContigOrdering(ref);

        traversalEngine = new TraverseLociByReference( reads, refFile, new java.util.ArrayList() );
    }


    public void execute( Walker walker,                        // the analysis technique to use.
                         List<GenomeLoc> locations ) {         // list of work to do
        ShardStrategy shardStrategy = null;
        if( locations != null )
            shardStrategy = ShardStrategyFactory.shatter( ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,
                                                          ref.getSequenceDictionary(),
                                                          SHARD_SIZE,
                                                          locations );
        else
            shardStrategy = ShardStrategyFactory.shatter( ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,
                                                          ref.getSequenceDictionary(),
                                                          SHARD_SIZE );

        ReferenceIterator refIter = new ReferenceIterator(ref);
        SAMDataSource dataSource = null;

        try {
            dataSource = new SAMDataSource( Arrays.asList( new String[] { reads.getCanonicalPath() } ) );
        }
        catch( SimpleDataSourceLoadException ex ) {
            throw new RuntimeException( ex );
        }        
        catch( IOException ex ) {
            throw new RuntimeException( ex );
        }

        Object accumulator = ((LocusWalker<?,?>)walker).reduceInit();

        for(Shard shard: shardStrategy) {
            Iterator<SAMRecord> readShard = null;
            try {
                readShard = dataSource.seek( shard.getGenomeLoc() );
            }
            catch( SimpleDataSourceLoadException ex ) {
                throw new RuntimeException( ex );
            }

            ReferenceProvider referenceProvider = new ReferenceProvider( refIter );
            LocusContextProvider locusProvider = new LocusContextProvider( readShard );

            accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );
        }

        traversalEngine.printOnTraversalDone("loci", accumulator);        
        walker.onTraversalDone(accumulator);        
    }

    
}
