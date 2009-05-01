package org.broadinstitute.sting.gatk.executive;

import java.io.OutputStream;
/**
 * User: hanna
 * Date: May 1, 2009
 * Time: 11:47:44 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A queueable task to merge two output files.  Provides a static 'lock' to
 * track whether a merge is currently happening.
 */
public class OutputMerger implements Runnable {
    public final ShardOutput shardOutput;

    private final OutputStream out;
    private final OutputStream err;

    /**
     * Indicates whether a merge has been queued.  Keeps merges from stepping on each other.
     */
    private static boolean mergeQueued;

    public OutputMerger( ShardOutput shardOutput, OutputStream out, OutputStream err ) {
        this.shardOutput = shardOutput;
        this.out = out;
        this.err = err;
    }

    public boolean isComplete() {
        return shardOutput.isComplete();
    }

    /**
     * Is a merge currently queued?
     * @return True if a merge is waiting on the queue.  False otherwise.
     */
    public static boolean isMergeQueued() {
        return mergeQueued;
    }

    public static void queueMerge() {
        if( mergeQueued )
            throw new IllegalStateException( "Attempted to mark that a merge has been queued when a merge is already queued." );
        mergeQueued = true;
    }

    public void run() {
        if( !mergeQueued )
            throw new IllegalStateException( "Starting a file merge, but our state shows no merge has been queued." );
        try {
            shardOutput.mergeInto( out, err );
        }
        finally {
            mergeQueued = false;
        }
    }
}
