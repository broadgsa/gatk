package org.broadinstitute.sting.gatk.walkers;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 26, 2009
 * Time: 5:34:11 PM
 * To change this template use File | Settings | File Templates.
 */

/**
 * Indicates that a class is tree reducible, aka that any two adjacent
 * shards of the data can reduce with each other, and the composite result
 * can be reduced with other composite results.
 */
public interface TreeReducible<ReduceType> {
    /**
     * A composite, 'reduce of reduces' function.
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    ReduceType treeReduce(ReduceType lhs, ReduceType rhs);
}
