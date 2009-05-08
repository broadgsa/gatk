package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.dataSources.providers.LocusContextProvider;
import org.broadinstitute.sting.gatk.dataSources.providers.ReferenceProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.GenomeAnalysisTK;
import org.broadinstitute.sting.gatk.OutputTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.concurrent.Callable;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
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
    private OutputMerger output;

    public ShardTraverser( TraverseLociByReference traversalEngine,
                           Walker walker,
                           Shard shard,
                           IndexedFastaSequenceFile reference,
                           SAMDataSource reads,
                           OutputMerger output ) {
        this.walker = walker;
        this.traversalEngine = traversalEngine;
        this.shard = shard;
        this.reference = reference;
        this.reads = reads;
        this.output = output;
    }

    public Object call() {
        Object accumulator = walker.reduceInit();

        CloseableIterator<SAMRecord> readShard = null;
        readShard = reads.seek( shard );

        ReferenceProvider referenceProvider = new ReferenceProvider( reference, shard );
        LocusContextProvider locusProvider = new LocusContextProvider( readShard );

        OutputTracker outputTracker = GenomeAnalysisTK.Instance.getOutputTracker();

        outputTracker.setLocalStreams( output.getOutStream(), output.getErrStream() );
        try {
            accumulator = traversalEngine.traverse( walker, shard, referenceProvider, locusProvider, accumulator );
        }
        finally {
            readShard.close();

            output.complete();
            outputTracker.removeLocalStreams();            
        }

        return accumulator;
    }
}
