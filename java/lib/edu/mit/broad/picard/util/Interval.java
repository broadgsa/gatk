package edu.mit.broad.picard.util;

import edu.mit.broad.picard.PicardException;

import java.util.List;
import java.util.Collection;

/**
 * Represents a simple interval on a sequence.  Coordinates are 1-based closed ended.
 *
 * @author Tim Fennell
 */
public class Interval implements Comparable<Interval>, Cloneable {
    private String sequence;
    private int start;
    private int end;
    private boolean negativeStrand;
    private String name;

    /**
     * Constructs an interval with the supplied sequence and start and end. If the end
     * position is less than the start position an exception is thrown.
     *
     * @param sequence the name of the sequence
     * @param start the start position of the interval on the sequence
     * @param end the end position of the interval on the sequence
     */
    public Interval(String sequence, int start, int end) {
        this.sequence = sequence;
        this.start = start;
        this.end = end;

        if (this.end < this.start) throw new IllegalArgumentException("start must be less than or equal to end!");
    }

    /**
     * Constructs an interval with the supplied sequence and start, end, strand and name.
     * If the end position is less than the start position an exception is thrown.
     *
     * @param sequence the name of the sequence
     * @param start the start position of the interval on the sequence
     * @param end the end position of the interval on the sequence
     * @param negative true to indicate negative strand, false otherwise
     * @param name the name (possibly null) of the interval
     *
     */
    public Interval(String sequence, int start, int end, boolean negative, String name) {
        this(sequence, start, end);
        this.negativeStrand = negative;
        this.name = name;
    }

    /** Gets the name of the sequence on which the interval resides. */
    public String getSequence() { return sequence; }

    /** Gets the 1-based start position of the interval on the sequence. */
    public int getStart() { return start; }

    /** Gets the 1-based closed-ended end position of the interval on the sequence. */
    public int getEnd() { return end; }

    /** Returns true if the interval is on the negative strand, otherwise false. */
    public boolean isNegativeStrand() { return this.negativeStrand; }

    /** Returns true if the interval is on the positive strand, otherwise false. */
    public boolean isPositiveStrand() { return !this.negativeStrand; }

    /** Returns the name of the interval, possibly null. */
    public String getName() { return this.name; }

    /** Returns true if this interval overlaps the other interval, otherwise false. */
    public boolean intersects(Interval other) {
        return  (this.getSequence().equals(other.getSequence()) &&
                 CoordMath.overlaps(this.start, this.end, other.start, other.end));
    }

    /** Returns true if this interval overlaps the other interval, otherwise false. */
    public boolean abuts(Interval other) {
        return this.getSequence().equals(other.getSequence()) &&
               (this.start == other.end + 1 || other.start == this.end + 1);
    }

    /** Gets the length of this interval. */
    public int length() { return this.end - this.start + 1; }

    /** Counts the total number of bases a collection of intervals. */
    public static long countBases(Collection<Interval> intervals) {
        long total = 0;
        for (Interval i : intervals) {
            total += i.length();
        }

        return total;
    }


    /**
     * Sort based on sequence.compareTo, then start pos, then end pos
     * with null objects coming lexically last
     */
    public int compareTo(Interval that) {
        if (that == null) return -1; // nulls last

        int result = this.getSequence().compareTo(that.getSequence());
        if (result == 0) {
            if (this.start == that.start) {
                result = this.end - that.end;
            }
            else {
                result = this.start - that.start;
            }
        }

        return result;
    }

    /** Equals method that agrees with {@link #compareTo(Interval)}. */
    public boolean equals(Interval that) {
        return (this.compareTo(that) == 0);
    }

    public int hashCode() {
        int result;
        result = sequence.hashCode();
        result = 31 * result + (start ^ (start >>> 32));
        result = 31 * result + (end ^ (end >>> 32));
        return result;
    }

    public String toString() {
        return getSequence() + ":" + start + "-" + end;
    }

    @Override
    public Interval clone() {
        try { return (Interval) super.clone(); }
        catch (CloneNotSupportedException cnse) { throw new PicardException("That's unpossible", cnse); }
    }
}
