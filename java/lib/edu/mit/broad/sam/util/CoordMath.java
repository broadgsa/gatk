/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.util;

public class CoordMath {

    public static long getLength(final long start, final long end) {
        return (end - start) + 1;
    }

    public static long getStart(final long end, final long length) {
        return end - length + 1;
    }

    public static long getEnd(final long start, final long length) {
        return start + length - 1;
    }

    /**
     * Offsets are meant to exclude the 'offset' number of bases
     */
    public static long getStartFromOffset(final long offset, final long length) {
        return offset + 1;
    }

    public static long getEndFromOffset(final long offset, final long length) {
        return length - offset;
    }

    public static long getLengthFromOffsets(final long startOffset, final long endOffset, final long length) {
        return getLength(getStartFromOffset(startOffset, length),
                         getEndFromOffset(endOffset, length));
    }

    /**
     * Gets a sub-sequence from a java.lang.String (which is zero based) using one based
     * sequence coordinated.  The base at the end coordinate will be included.
     *
     * @param sequence The String of base pairs
     * @param begin The one based start coordinate
     * @param end The one based end coordinate
     * @return The subsequence specified
     */
    public static String getSubsequence(final String sequence, final int begin, final int end) {
        return sequence.substring(begin-1, end);
    }

    /**
     * Checks to see if the two sets of coordinates have any overlap.
     */
    public static boolean overlaps(final long start, final long end, final long start2, final long end2) {
        return (start2 >= start && start2 <= end) || (end2 >=start && end2 <= end) ||
                encloses(start2, end2, start, end);
    }

    /** Returns true if the "inner" coords and totally enclosed by the "outer" coords. */
    public static boolean encloses(final long outerStart, final long outerEnd, final long innerStart, final long innerEnd) {
        return innerStart >= outerStart && innerEnd <= outerEnd;
    }

    /**
     * Determines the amount of overlap between two coordinate ranges. Assumes that the two ranges
     * actually do overlap and therefore may produce strange results when they do not!
     */
    public static long getOverlap(final long start, final long end, final long start2, final long end2) {
        return getLength(Math.max(start, start2), Math.min(end, end2));
    }
}
