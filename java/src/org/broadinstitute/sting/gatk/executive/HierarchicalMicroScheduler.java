package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.traversals.TraverseLociByReference;
import org.broadinstitute.sting.gatk.traversals.TraversalEngine;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 5:41:04 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * A microscheduler that schedules shards according to a tree-like structure.
 * Requires a special walker tagged with a 'TreeReducible' interface.
 */
public class HierarchicalMicroScheduler extends MicroScheduler {
    /**
     * How many threads should the hierarchical scheduler try to keep busy.
     */
    private int nThreadsToUse;

    private TraverseLociByReference traversalEngine = null;

    /**
     * Create a new hierarchical microscheduler to process the given reads and reference.
     * @param reads Reads file(s) to process.
     * @param refFile Reference for driving the traversal.
     * @param nThreadsToUse maximum number of threads to use to do the work
     */
    protected HierarchicalMicroScheduler( List<File> reads, File refFile, int nThreadsToUse ) {
        super( reads, refFile );
        this.nThreadsToUse = nThreadsToUse;
        traversalEngine = new TraverseLociByReference( reads, refFile, new java.util.ArrayList() );
    }

    public TraversalEngine getTraversalEngine() {
        return traversalEngine;
    }

    public void execute( Walker walker, List<GenomeLoc> intervals ) {
        
    }
}
