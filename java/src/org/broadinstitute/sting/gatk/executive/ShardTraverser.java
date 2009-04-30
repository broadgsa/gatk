package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.concurrent.Callable;
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 4:40:38 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */
/**
 * Carries the walker over a given shard, in a callable interface.
 */
public class ShardTraverser implements Callable {
    private Walker walker;
    private TraverseLociByReference traversalEngine;
    private Shard shard;
    private IndexedFastaSequenceFile reference;
    private SAMDataSource reads;

    public ShardTraverser( TraverseLociByReference traversalEngine,
                           Walker walker,
                           Shard shard,
                           IndexedFastaSequenceFile reference,
                           SAMDataSource reads ) {
        this.walker = walker;
        this.traversalEngine = traversalEngine;
        this.shard = shard;
        this.reference = reference;
        this.reads = reads;
    }

    public Object call() {
        Object accumulator = ((LocusWalker<?,?>)walker).reduceInit();

        MergingSamRecordIterator2 readShard = null;
        try {
            readShard = (MergingSamRecordIterator2)reads.seek( shard );
        }
        catch( SimpleDataSourceLoadException ex ) {
            throw new RuntimeException( ex );
        }

        ReferenceProvider referenceProvider = new ReferenceProvider( reference, shard.getGenomeLoc() );
        LocusContextProvider locusProvider = new LocusContextProvider( readShard );

        accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );

        readShard.close();

        return accumulator;
    }
}
