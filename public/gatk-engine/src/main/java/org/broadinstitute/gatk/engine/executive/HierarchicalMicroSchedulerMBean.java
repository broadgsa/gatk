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
/**
 * User: hanna
 * Date: May 29, 2009
 * Time: 4:05:27 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * An interface for retrieving runtime statistics about how the hierarchical
 * microscheduler is behaving. 
 */
public interface HierarchicalMicroSchedulerMBean {
    /**
     * How many tree reduces are waiting in the tree reduce queue?
     * @return Total number of reduces waiting in the tree reduce queue?
     */
    public int getNumberOfTasksInReduceQueue();

    /**
     * How many pending I/O combining tasks are waiting in the queue?
     * @return Total number of I/O tasks waiting in the I/O queue.
     */
    public int getNumberOfTasksInIOQueue();

    /**
     * What is the total time spent running traversals?
     * @return Total time spent traversing shards; 0 if none have been traversed.
     */
    public long getTotalShardTraverseTimeMillis();

    /**
     * What is the average time spent running traversals?
     * @return Average time spent traversing shards; 0 if none have been traversed.
     */
    public long getAvgShardTraverseTimeMillis();

    /**
     * What is the total time spent merging output?
     */
    public long getTotalOutputMergeTimeMillis();

    /**
     * What is the total time spent running tree reduces?
     * @return Total time spent running tree reduces; 0 if none have been run.
     */
    public long getTotalTreeReduceTimeMillis();

    /**
     * What is the average time spent running tree reduces?
     * @return Average time spent running tree reduces; 0 if none have been run.
     */
    public long getAvgTreeReduceTimeMillis();
}

