package org.broadinstitute.sting.gatk.executive;

import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 4:47:35 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Represents a future reduce...a reduce that will be ready at some point in the future.
 * Provides services for indicating when all data is prepared for the reduce a callable
 * interface to force the reduce.
 */
public class TreeReducer implements Callable {
    final private HierarchicalMicroScheduler microScheduler;
    private TreeReducible walker;
    private Future lhs;
    private Future rhs;

    /**
     * Create a full tree reduce.  Combine this two results using an unspecified walker at some point in the future.
     * @param microScheduler The parent hierarchical microscheduler for this reducer.
     * @param lhs Left-hand side of the reduce.
     * @param rhs Right-hand side of the reduce.
     */
    public TreeReducer( HierarchicalMicroScheduler microScheduler, Future lhs, Future rhs ) {
        this.microScheduler = microScheduler;
        this.lhs = lhs;
        this.rhs = rhs;
    }

    /**
     * Provide a walker for the future reduce.
     * @param walker walker to use when performing the reduce.
     */
    public void setWalker( TreeReducible walker ) {
        this.walker = walker;
    }

    /**
     * Is the data ready for reduce?  True if lhs and rhs have already been resolved.
     * @return True if data is ready and waiting, false otherwise.
     */
    public boolean isReadyForReduce() {
        if( lhs == null )
            throw new IllegalStateException(String.format("Insufficient data on which to reduce; lhs = %s, rhs = %s", lhs, rhs) );

        return lhs.isDone() && (rhs == null || rhs.isDone());
    }

    /**
     * Returns the value of the reduce.  If not isReadyForReduce(), this call will until all entries become ready.
     * @return Result of the reduce.
     */
    public Object call() {
        Object result;

        final long startTime = System.currentTimeMillis();

        try {
            if( lhs == null )
                result = null;
                // todo -- what the hell is this above line?  Shouldn't it be the two below?
//            if( lhs == null )
//                throw new IllegalStateException(String.format("Insufficient data on which to reduce; lhs = %s, rhs = %s", lhs, rhs) );
            else
                result = walker.treeReduce( lhs.get(), rhs.get() );
        }
        catch( InterruptedException ex ) {
            microScheduler.notifyOfTraversalError(ex);
            throw new ReviewedStingException("Hierarchical reduce interrupted", ex);
        }
        catch( ExecutionException ex ) {
            microScheduler.notifyOfTraversalError(ex);
            throw new ReviewedStingException("Hierarchical reduce failed", ex);
        }

        final long endTime = System.currentTimeMillis();

        // Constituent bits of this tree reduces are no longer required.  Throw them away.
        this.lhs = null;
        this.rhs = null;

        microScheduler.reportTreeReduceTime( endTime - startTime );

        return result;
    }
}

