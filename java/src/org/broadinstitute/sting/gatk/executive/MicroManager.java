package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * A micro-scheduling manager for N-way threaded execution of a traversal
 *
 */
public class MicroManager {

    public MicroManager( TraversalEngineExecutive TEfactory,  // makes worker units
                         GenomeLoc[] locations,             // list of work to do
                         int nThreadsToUse,                 // maximum number of threads to use to do the work
                         int initialChunkSize ) {           // the initial chunk size for dividing up the work
        // do a lot of work here to organize the computation
    }

    public void execute() {
        // actually divide up the work, create worker units, and send them off to do work

    }
}
