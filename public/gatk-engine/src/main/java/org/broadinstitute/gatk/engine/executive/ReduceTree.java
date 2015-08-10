/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.executive;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.Future;
/**
 * User: hanna
 * Date: Apr 28, 2009
 * Time: 11:09:29 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A tree for organizing reduce results and detecting when enough dependencies
 * are resolved for a reduce to be scheduled.  The tree can trigger a callback
 * whenever it believes a reduce operation is pending.
 *
 * Not thread-safe.  All calls should be made sequentially from the same thread.
 */
public class ReduceTree {
    /**
     * Data structure for the tree.  Each entry in the outer list represents a level
     * of the tree, and each entry in the inner queues represent nodes in that level.
     *
     * Whenever a reduce can happen, the entries to be reduced are pulled out of
     * their slots in level n of the tree and the composite entry is added to level n+1.
     */
    private List<Queue<Future>> treeNodes = new ArrayList<Queue<Future>>();

    /**
     * The entire contents have been added to the tree.  Completely schedule the reductions.
     */
    private boolean treeComplete = false;

    /**
     * Called to indicate that all data required to perform a given reduce has been scheduled.
     */
    private TreeReduceNotifier treeReduceNotifier = null;

    /**
     * Creates a ReduceTree.
     * @param notifier A callback indicating that all data required to perform a given reduce has been scheduled.
     */
    public ReduceTree( TreeReduceNotifier notifier ) {
        this.treeReduceNotifier = notifier;
    }

    /**
     * A callback indicating that all computations have been scheduled to complete the given reduce.
     */
    public interface TreeReduceNotifier {
        /**
         * Indicates that a reduce is ready to happen.
         * @param lhs Left-hand side of the tree reduce.
         * @param rhs Right-hand side of the tree reduce.
         * @return The future result of the computation reduce(lhs,rhs)
         */
        Future notifyReduce( Future lhs, Future rhs );
    }

    /**
     * Add an entry to the list of data to be reduced.  The results of entry.get() will
     * be scheduled for reduction with neighboring elements.
     * @param entry Entry to be paired with other elements.
     */
    public void addEntry( Future entry ) {
        addNodeAtLevel( entry, 0 );
    }

    /**
     * Signal to the ReduceTree that all possible data has been added and it should reduce
     * as much as is possible.
     */
    public void complete() {
        treeComplete = true;
        reduce();
    }

    /**
     * Gets the placeholder for the final result of the tree reduce.
     * @return Future whose get() method will return the final result.  Null if nothing has been added.
     */
    public Future getResult() {
        if( !treeComplete )
            throw new IllegalStateException( "Cannot get the final result for an incomplete tree.");

        // If nothing has been added to the tree, return null.
        if( treeNodes.size() == 0 )
            return null;

        // Assert that there aren't any pending computations that were forgotten along the way.
        for( int i = 0; i < treeNodes.size() - 2; i++ ) {
            if( treeNodes.get(i).size() > 0 )
                throw new IllegalStateException( "Some inner reduces were missed along the way.");
        }

        Queue<Future> lastLevel = treeNodes.get(treeNodes.size() - 1);

        // Assert that there's only one reduce left at the last level.
        if( lastLevel.size() != 1 )
            throw new IllegalStateException( "Invalid number of entries at the tip of the tree: " + lastLevel.size() );

        // Get the placeholder for the last result.
        return lastLevel.element();
    }

    /**
     * Recursively collapse the tree whereever possible.
     */
    protected void reduce() {
        reduce( 0 );
    }

    /**
     * Recursively collapse the tree, starting at the specified level.
     * @param level Level at which to start reducing.
     */
    private void reduce( int level ) {
        // base case for recursion.
        if( treeNodes.size() <= level )
            return;

        Queue<Future> treeLevel = treeNodes.get(level);

        while( treeLevel.size() >= 2 ) {
            addNodeAtLevel( treeReduceNotifier.notifyReduce( treeLevel.remove(), treeLevel.remove() ), level + 1 );
        }

        if( treeLevel.size() == 1 && treeComplete && !isDeepestLevel(level) ) {
            Future element = treeLevel.remove();
            addNodeAtLevel( element, level + 1 );
        }

        reduce( level + 1 );
    }

    private boolean isDeepestLevel( int level ) {
        return level == (treeNodes.size() - 1);
    }

    /**
     * Add the given node to the tree at the corresponding level.  Create the level
     * if it doesn't exist.
     * @param node Node to add.  Must not be null.
     * @param level Level number at which to add.  0-based index into treeNodes list.
     */
    protected void addNodeAtLevel( Future node, int level ) {
        while( treeNodes.size() <= level )
            treeNodes.add( new LinkedList<Future>() );
        treeNodes.get(level).add(node);
        reduce(level);
    }

}
