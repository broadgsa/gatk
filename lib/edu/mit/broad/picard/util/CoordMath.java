/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright 2005 by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;


/**
 * Basic coordinate-based math utils, so it's encapsulated in one place!  Assumes
 * a one-based coordinate system and then 'end' is always inclusive
 */
public class CoordMath {

    /** Gets the length of an interval given the start and the end. */
    public static int getLength(int start, int end) { return (end - start) + 1; }

    /** Gets the start of an interval given the end and the length. */
    public static int getStart(int end, int length) { return end - length + 1; }

    /** Gets the end of an interval given the start and the length. */
    public static int getEnd(int start, int length) { return start + length - 1; }

    /** Checks to see if the two sets of coordinates have any overlap. */
    public static boolean overlaps(int start, int end, int start2, int end2) {
        return (start2 >= start && start2 <= end) || (end2 >=start && end2 <= end) ||
                encloses(start2, end2, start, end);
    }

    /** Returns true if the "inner" coords and totally enclosed by the "outer" coords. */
    public static boolean encloses(int outerStart, int outerEnd, int innerStart, int innerEnd) {
        return innerStart >= outerStart && innerEnd <= outerEnd;
    }

    /**
     * Determines the amount of overlap between two coordinate ranges. Assumes that the two ranges
     * actually do overlap and therefore may produce strange results when they do not!
     */
    public static int getOverlap(int start, int end, int start2, int end2) {
        return getLength(Math.max(start, start2), Math.min(end, end2));
    }

    /** 
     * Determines the read cycle number for the base
     * 
     *  @param isNegativeStrand true if the read is negative strand
     *  @param readLength
     *  @param readBaseIndex the 0-based index of the read base in question
     */
    public static int getCycle(boolean isNegativeStrand, int readLength, final int readBaseIndex) {
        return isNegativeStrand ? readLength - readBaseIndex : readBaseIndex + 1;
    }
}
