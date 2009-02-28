/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.util;

public class SequenceUtil {
    /**
     * Calculate the reverse complement of the specified sequence
     * (Stolen from Reseq)
     *
     * @param sequenceData
     * @return reverse complement
     */
    public static String reverseComplement(String sequenceData) {

        final char[] original = sequenceData.toCharArray();
        final char[] complement = new char[original.length];

        for (int i=0, j=complement.length-1; i<original.length; ++i, --j) {
            switch ( original[i] ) {

                // inlined for performance (although HotSpot may do this anyway...)
                case 'a' : complement[j] = 't'; break;
                case 'A' : complement[j] = 'T'; break;
                case 'c' : complement[j] = 'g'; break;
                case 'C' : complement[j] = 'G'; break;
                case 'g' : complement[j] = 'c'; break;
                case 'G' : complement[j] = 'C'; break;
                case 't' : complement[j] = 'a'; break;
                case 'T' : complement[j] = 'A'; break;
                default  : complement[j] = original[i];
            }
        }

        return new String(complement);
    }

    /** Attempts to efficiently compare two bases stored as bytes for equality. */
    public static boolean basesEqual(byte lhs, byte rhs) {
        if (lhs == rhs) return true;
        else {
            if (lhs > 90) lhs -= 32;
            if (rhs > 90) rhs -= 32;
        }

        return lhs == rhs;
    }
    
    /**
     * returns true if the value of base represents a no call
     */
    public static boolean isNoCall(byte base) {
        return base == 'N' || base == 'n' || base == '.';
    }

}
