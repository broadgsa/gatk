package org.broadinstitute.sting.utils;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 19, 2009
 * Time: 12:03:39 PM
 * To change this template use File | Settings | File Templates.
 */

/** Abstraction of a closed interval [start,stop] 
 *
 */
public interface Interval {
    /** Start position of the interval.
     *
     * @return <start> for the interval [start,stop]
     */
    public long getStart();

    /** Sets start position of the interval.
     *
     * @param s start coordinate
     */
    public void setStart(long s);

    /** End position of the interval.
     *
     * @return <stop> for the interval [start,stop]
     */
    public long getStop();

    /** Sets stop position of the interval.
     *
     * @param s stop coordinate
     */
    public void setStop(long s);

    /** Length of the interval. There is currently no contract, an implementation may return negative length
     * or a length inconsistent with getStop() - getStart() + 1 if it chooses so.
     *
     * @return a number representing the length of the interval according to specific implementation
     */
    public long getLength();

    /** Returns true if this interval overlaps with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
      * @param i Another interval
     * @return true iff intervals overlap
     */
    public boolean overlapsP(org.broadinstitute.sting.utils.Interval i);

    /** Returns true if this interval does not overlap with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
      * @param i Another interval
     * @return true iff intervals do not overlap
     */
    public boolean disjointP(org.broadinstitute.sting.utils.Interval i);
}
