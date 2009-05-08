package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.dataSources.providers.ShardDataProvider;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.GenomeAnalysisTK;
import org.broadinstitute.sting.gatk.OutputTracker;
import org.broadinstitute.sting.gatk.walkers.Walker;

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
    private ShardDataProvider dataProvider;
    private OutputMerger output;

    public ShardTraverser( TraverseLociByReference traversalEngine,
                           Walker walker,
                           Shard shard,
                           ShardDataProvider dataProvider,
                           OutputMerger output ) {
        this.walker = walker;
        this.traversalEngine = traversalEngine;
        this.shard = shard;
        this.dataProvider = dataProvider;
        this.output = output;
    }

    public Object call() {
        Object accumulator = walker.reduceInit();
        OutputTracker outputTracker = GenomeAnalysisTK.Instance.getOutputTracker();
        outputTracker.setLocalStreams( output.getOutStream(), output.getErrStream() );

        try {
            accumulator = traversalEngine.traverse( walker, shard, dataProvider, accumulator );
        }
        finally {
            dataProvider.close();
            output.complete();
            outputTracker.removeLocalStreams();            
        }

        return accumulator;
    }
}
