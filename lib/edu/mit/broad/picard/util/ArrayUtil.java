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

public class ArrayUtil {

    /**
     * Reverse the elements of the given array in place
     */
    public static <T> void reverseArray(T[] array) {
        for (int left=0, right=array.length-1; left<right; left++, right--) {
            // exchange the first and last
            T temp = array[left]; array[left]  = array[right]; array[right] = temp;
        }
    }

    /**
     * clone the above method as necessary for non-object types
     */
    public static void reverseArray(byte[] array) {
        for (int left=0, right=array.length-1; left<right; left++, right--) {
            // exchange the first and last
            byte temp = array[left]; array[left]  = array[right]; array[right] = temp;
        }
    }
}
