package org.broadinstitute.sting.gatk.executive;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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

    protected List<GenomeLoc> intervalList = null;

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

    public void setIntervalList(List<GenomeLoc> intervalList) {
        this.intervalList = intervalList;
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
            // todo: remove this code when we acutally handle command line args of multiple bam files
            ArrayList<File> fl = new ArrayList<File>();
            if (reads.getName().endsWith(".list")) {
                BufferedReader bis = new BufferedReader(new FileReader(reads));
                String line = null;
                while ((line = bis.readLine()) != null) {
                    if (!line.equals("")){
                        fl.add(new File(line));
                    }
                }

            } else {
                fl.add(reads);
            }
            dataSource = new SAMDataSource( fl );
        }
        catch( SimpleDataSourceLoadException ex ) {
            throw new RuntimeException( ex );
        }
        catch( IOException ex ) {
            throw new RuntimeException( ex );
        }

        Object accumulator = ((LocusWalker<?,?>)walker).reduceInit();

        for(Shard shard: shardStrategy) {
            // CloseableIterator<SAMRecord> readShard = null;
            MergingSamRecordIterator2 readShard = null;
            try {
                readShard = dataSource.seek( shard.getGenomeLoc() );
            }
            catch( SimpleDataSourceLoadException ex ) {
                throw new RuntimeException( ex );
            }

            ReferenceProvider referenceProvider = new ReferenceProvider( refIter );
            LocusContextProvider locusProvider = new LocusContextProvider( readShard );

            accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );
            readShard.close();
        }

        traversalEngine.printOnTraversalDone("loci", accumulator);
        walker.onTraversalDone(accumulator);
    }


}
