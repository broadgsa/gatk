/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import edu.mit.broad.sam.util.CoordMath;

class SimpleMapping implements Comparable<SimpleMapping> {
    String arachneIndex;
    long startPos;
    long endPos;
    String sequenceName;

    public SimpleMapping(final String arachneIndex, final long startPos, final long endPos, final String sequenceName) {
        this.arachneIndex = arachneIndex;
        this.startPos = startPos;
        this.endPos = endPos;
        this.sequenceName = sequenceName;

        if (this.endPos < this.startPos) throw new IllegalArgumentException("startPos must be less than endPos!");
    }

    public String getArachneIndex() {
        return arachneIndex;
    }

    public void setArachneIndex(final String arachneIndex) {
        this.arachneIndex = arachneIndex;
    }

    public long getStartPos() {
        return startPos;
    }

    public void setStartPos(final long startPos) {
        this.startPos = startPos;
    }

    public long getEndPos() {
        return endPos;
    }

    public void setEndPos(final long endPos) {
        this.endPos = endPos;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public void setSequenceName(final String sequenceName) {
        this.sequenceName = sequenceName;
    }

    public SimpleMapping intersection(final SimpleMapping other) {
        if (this.intersects(other)) {
            return new SimpleMapping(this.getArachneIndex(),
                    (this.getStartPos() >= other.getStartPos())?this.getStartPos():other.getStartPos(),
                    (this.getEndPos() <= other.getEndPos())?this.getEndPos():other.getEndPos(), this.getSequenceName());
        }

        return null;
    }

    public boolean intersects(final SimpleMapping other) {
        return  (this.getArachneIndex().equals(other.getArachneIndex()) &&
                CoordMath.overlaps(this.getStartPos(), this.getEndPos(), other.getStartPos(), other.getEndPos()));
    }

    public long length() {
        return CoordMath.getLength(startPos, endPos);
    }

    /**
     * Sort based on sequence.compareTo, then start pos, then end pos
     * with null objects coming lexically last
     */
    public int compareTo(final SimpleMapping that) {
        if (that == null) return -1; // nulls last

        int result = this.getArachneIndex().compareTo(that.getArachneIndex());
        if (result == 0) {
            if (this.getStartPos() == that.getStartPos()) {
                result = ((int) (this.getEndPos() - that.getEndPos()));
            } else {
                result = ((int) (this.getStartPos() - that.getStartPos()));
            }
        }

        // normalize to -1, 0, 1
        if (result > 1) result = 1;
        else if (result < -1) result = -1;
        return result;
    }

    public boolean equals(final SimpleMapping that) {
        return (this.compareTo(that) == 0);
    }

    public int hashCode() {
        int result;
        result = arachneIndex.hashCode();
        result = 31 * result + (int) (startPos ^ (startPos >>> 32));
        result = 31 * result + (int) (endPos ^ (endPos >>> 32));
        return result;
    }

    public String toString() {
        return getArachneIndex() + ":" + getStartPos() + "-" + getEndPos();
    }
}
