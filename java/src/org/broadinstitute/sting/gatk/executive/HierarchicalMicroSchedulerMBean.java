package org.broadinstitute.sting.gatk.executive;
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
public interface HierarchicalMicroSchedulerMBean extends MicroSchedulerMBean {
    /**
     * What is the total number of shards assigned to this microscheduler?
     * @return Total number of shards to process.
     */
    public int getTotalNumberOfShards();

    /**
     * How many shards are remaining for this microscheduler to process?
     * @return Remaining number of shards to process.
     */
    public int getRemainingNumberOfShards();

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

