package org.broadinstitute.sting.indels;

import org.broadinstitute.sting.utils.Interval;

/**  This class represents an indel as an interval with respect to the <i>original</i> reference and, in addition,
 * stores the indel type ( (I)nsertion or (D)eletion ) and can return meaningful event size (see below).
 * Depending on the indel type, the positions on the reference are:
 *    <ul>
 *    <li> Deletion ( e.g. deletion of ACAC from the ref: ACGTT[ACAC]TTTAG to ACGTT[]TTTAG) - start and stop of
 *         the interval are first and last deleted bases on the original reference (those in square brackets in the
 *         first sequence in the example)
 *    <li> Insertion ( e.g. insertion of GTGT into the ref: ACGTT{}TTTAG to ACGTT{GTGT}TTTAG) - start is the first
 *         position on the original reference <i>after</i> the insertion site (after the '}'), and stop is the last
 *         position <i>before</i> the insertion site (prior to '{').
 *    </ul>
 *
 * Given these definitions, the length of the interval, as returned by getLength() has the meaning of the length of
 * the event (affected bases) on the original reference: number of deleted bases for deletion and zero for insertion.
 * The length of the indel itself is returned by getIndelLength(), which is equal to getLength() for deletions and to
 * the actual number of inserted bases for insertions (while length on the reference, as returned by getLength() is zero).
 *
 * The overlaps are also meaningful with the above definitions: if an alignment to (or, in general, an interval on)
 * the original reference ends prior to <code>start</code>, or starts after <code>stop</code>, it does not overlap
 * with the indel event (neither spans over deleted region or contains any of the inserted bases).
 *
 */
public class Indel implements Interval {

    public static enum IndelType { I, D };
	
	private long mStart;
    private long mLength;
    private IndelType mType;

    /** Creates nBases-long indel at specified start position; the object will be unusable
     * until indel type is set.
     * @param start start position on the reference
     * @param nBases number of inserted or deleted bases
     */
    //public Indel(long start, long nBases) {
   //     mType=null;
   //     mStart=start;
   //     mLength=nBases;
   // }

    /** Creates nBases-long indel of the specified type (insertion or deletion), at specified start position.
     * @param start start position on the reference
     * @param nBases number of inserted or deleted bases
     * @param type Indel type: I or D.
     */
    public Indel(long start, long nBases, IndelType type) {
        mType=type;
        mStart=start;
        mLength=nBases;
    }

	/** Start coordinate on the reference; for deletions it is the position of the first deleted base,
	 * for insertions it is the first base after the insertion.
	 * This is the "left boundary" of the event on the original reference: every alignment that ends
     * befor this position on the reference does not overlap with the indel.
	 * @return indel's left boundary
	 */
	public long getStart() { return mStart; }

    /** Sets start position of the interval.
     *
     * @param s start coordinate
     */
	public void setStart(long s) { mStart = s; }

	/** Indel's stop coordinate on the reference; for deletions it is the position of the last deleted base,
	 * for insertions it is the last base before the insertion site (which makes it equal to getStart() - 1).
	 * This is the "right boundary" of the event: every alignment that starts after
     * this position on the reference
	 * does not overlap with the indel.
	 * @return indel's right boundary
	 */
	public long getStop() {
        if ( mType == IndelType.I ) return mStart - 1;
		else return mStart + mLength - 1;
	}

    /** This method is not supported in IndelInterval and will throw an exception. Use setIndelLength() instead.
     *
     * @param s stop coordinate
     */
    public void setStop(long s) {
        throw new UnsupportedOperationException("Method setStop(long) is not supported in IndelInterval");
    }

    /** Returns type of this indel ( I or D).
     *
     * @return I or D enum element
     */
	public IndelType getType() { return mType; }

    /** Sets the number of bases in this indel (i.e. the actual number of inserted or
     * deleted bases). Stop position will be always correctly computed based on the indel length and indel type.
     * @param nBases length of the indel (<i>not</i> the length of the event on the original reference!)
     */
	public void setIndelLength(long nBases) { mLength = nBases; }

    /** Returns actual number of inserted or deleted bases in the indel.
     *
     * @return number of bases (<i>not</i> the event length on the original reference).
     * @see #getLength()
     */
	public long getIndelLength() { return mLength; }

    /**
     * Returns true if this interval overlaps with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
     *
     * @param i Another interval
     * @return true iff intervals overlap
     */
    @Override
    public boolean overlapsP(Interval i) {
        return ! disjointP(i);  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Returns true if this interval does not overlap with i as judjed by getStart() and getStop() positions of the
     * two interval objects.
     *
     * @param i Another interval
     * @return true iff intervals do not overlap
     */
    @Override
    public boolean disjointP(Interval i) {
        return i.getStop() < this.getStart() || i.getStart() > this.getStop();
    }

    /** Returns length of the region affected by the indel on the original reference. Note that an insertion
	 *  has length of 0.
     *  @return length of the event on the original, unmodified reference
	 */
    @Override
	public long getLength() {
		if ( mType == IndelType.I ) return 0; 
		return mLength;
	}

    @Override
    public boolean equals(Object o) {
        if ( ! ( o instanceof Indel ) ) return false;
        Indel i = (Indel)o;
        return this.mType == i.mType && this.mStart == i.mStart && this.mLength == i.mLength ;
    }

    @Override
    public int hashCode() {
        return (int)( mStart << 2 + mLength );
    }
}
