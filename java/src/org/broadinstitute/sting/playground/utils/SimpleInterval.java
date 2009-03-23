package org.broadinstitute.sting.playground.utils;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 19, 2009
 * Time: 12:44:37 PM
 * To change this template use File | Settings | File Templates.
 */

/** Provides minimum complete implementation of Interval interface.
 */
public class SimpleInterval implements Interval {

    private long m_start;
    private long m_stop;

    public SimpleInterval(long start, long stop) { m_start = start; m_stop = stop; }

    public SimpleInterval(Interval i) { m_start = i.getStart(); m_stop = i.getStop(); }
//    public SimpleInterval() { m_start = -1; m_stop = -2; }

    /** Start position of the interval.
     *
     * @return <start> for the interval [start,stop]
     */
    @Override
    public long getStart() { return m_start; }

    /** Sets start position of the interval.
     *
     * @param s start coordinate
     */
    public void setStart(long s) { m_start = s; }

    /** End position of the interval.
     *
     * @return <stop> for the interval [start,stop]
     */
    @Override
    public long getStop() { return m_stop; }

    /** Sets stop position of the interval.
     *
     * @param s stop coordinate
     */
    public void setStop(long s) { m_stop = s; }

    /** Length of the interval. This default implementation returns getStop() - getStart() + 1.
     *
     * @return
     */
    @Override
    public long getLength() { return (m_stop - m_start + 1); };

    /** Returns true if this interval overlaps with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
      * @param i Another interval
     * @return true iff intervals overlap
     */
    @Override
    public boolean overlapsP(org.broadinstitute.sting.playground.utils.Interval i) {
        return ! disjointP(i);
    }

    /** Returns true if this interval does not overlap with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
      * @param i Another interval
     * @return true iff intervals do not overlap
     */
    @Override
    public boolean disjointP(org.broadinstitute.sting.playground.utils.Interval i) {
        return ( i.getStop() < this.m_start || i.getStart() > this.m_stop );
    }
}
