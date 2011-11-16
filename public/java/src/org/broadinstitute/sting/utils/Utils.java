/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import net.sf.samtools.util.StringUtil;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.Pair;

import java.net.InetAddress;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:12:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class Utils {
    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(Utils.class);

    public static final float JAVA_DEFAULT_HASH_LOAD_FACTOR = 0.75f;

    /**
     * Calculates the optimum initial size for a hash table given the maximum number
     * of elements it will need to hold. The optimum size is the smallest size that
     * is guaranteed not to result in any rehash/table-resize operations.
     *
     * @param maxElements  The maximum number of elements you expect the hash table
     *                     will need to hold
     * @return             The optimum initial size for the table, given maxElements
     */
    public static int optimumHashSize ( int maxElements ) {
        return (int)(maxElements / JAVA_DEFAULT_HASH_LOAD_FACTOR) + 2;
    }

    /**
     * Compares two objects, either of which might be null.
     *
     * @param lhs One object to compare.
     * @param rhs The other object to compare.
     *
     * @return True if the two objects are equal, false otherwise.
     */
    public static boolean equals(Object lhs, Object rhs) {
        if (lhs == null && rhs == null) return true;
        else if (lhs == null) return false;
        else return lhs.equals(rhs);
    }

    public static <T> List<T> cons(final T elt, final List<T> l) {
        List<T> l2 = new ArrayList<T>();
        l2.add(elt);
        if (l != null) l2.addAll(l);
        return l2;
    }

    public static void warnUser(final String msg) {
        warnUser(logger, msg);
    }
    
    public static void warnUser(final Logger logger, final String msg) {
        logger.warn(String.format("********************************************************************************"));
        logger.warn(String.format("* WARNING:"));
        logger.warn(String.format("*"));
        prettyPrintWarningMessage(logger, msg);
        logger.warn(String.format("********************************************************************************"));
    }

    /**
     * pretty print the warning message supplied
     *
     * @param logger logger for the message
     * @param message the message
     */
    private static void prettyPrintWarningMessage(Logger logger, String message) {
        StringBuilder builder = new StringBuilder(message);
        while (builder.length() > 70) {
            int space = builder.lastIndexOf(" ", 70);
            if (space <= 0) space = 70;
            logger.warn(String.format("* %s", builder.substring(0, space)));
            builder.delete(0, space + 1);
        }
        logger.warn(String.format("* %s", builder));
    }

    public static ArrayList<Byte> subseq(char[] fullArray) {
        byte[] fullByteArray = new byte[fullArray.length];
        StringUtil.charsToBytes(fullArray, 0, fullArray.length, fullByteArray, 0);
        return subseq(fullByteArray);
    }

    public static ArrayList<Byte> subseq(byte[] fullArray) {
        return subseq(fullArray, 0, fullArray.length - 1);
    }

    public static ArrayList<Byte> subseq(byte[] fullArray, int start, int end) {
        assert end < fullArray.length;
        ArrayList<Byte> dest = new ArrayList<Byte>(end - start + 1);
        for (int i = start; i <= end; i++) {
            dest.add(fullArray[i]);
        }
        return dest;
    }

    public static String baseList2string(List<Byte> bases) {
        byte[] basesAsbytes = new byte[bases.size()];
        int i = 0;
        for (Byte b : bases) {
            basesAsbytes[i] = b;
            i++;
        }
        return new String(basesAsbytes);
    }

    /**
     * join the key value pairs of a map into one string, i.e. myMap = [A->1,B->2,C->3] with a call of:
     * joinMap("-","*",myMap) -> returns A-1*B-2*C-3
     *
     * Be forewarned, if you're not using a map that is aware of the ordering (i.e. HashMap instead of LinkedHashMap)
     * the ordering of the string you get back might not be what you expect! (i.e. C-3*A-1*B-2 vrs A-1*B-2*C-3)
     *
     * @param keyValueSeperator the string to seperate the key-value pairs
     * @param recordSeperator the string to use to seperate each key-value pair from other key-value pairs
     * @param map the map to draw from
     * @param <L> the map's key type
     * @param <R> the map's value type
     * @return a string representing the joined map
     */
    public static <L,R> String joinMap(String keyValueSeperator, String recordSeperator, Map<L,R> map) {
        if (map.size() < 1) { return null; }
        String joinedKeyValues[] = new String[map.size()];
        int index = 0;
        for (L key : map.keySet()) {
           joinedKeyValues[index++] = String.format("%s%s%s",key.toString(),keyValueSeperator,map.get(key).toString());
        }
        return join(recordSeperator,joinedKeyValues);
    }

    /**
     * Splits a String using indexOf instead of regex to speed things up.
     *
     * @param str the string to split.
     * @param delimiter the delimiter used to split the string.
     * @return an array of tokens.
     */
    public static ArrayList<String> split(String str, String delimiter) {
        return split(str, delimiter, 10);
    }

    /**
     * Splits a String using indexOf instead of regex to speed things up.
     *
     * @param str the string to split.
     * @param delimiter the delimiter used to split the string.
     * @param expectedNumTokens The number of tokens expected. This is used to initialize the ArrayList.
     * @return an array of tokens.
     */
    public static ArrayList<String> split(String str, String delimiter, int expectedNumTokens) {
        final ArrayList<String> result =  new ArrayList<String>(expectedNumTokens);

        int delimiterIdx = -1;
        do {
            final int tokenStartIdx = delimiterIdx + 1;
            delimiterIdx = str.indexOf(delimiter, tokenStartIdx);
            final String token = (delimiterIdx != -1 ? str.substring(tokenStartIdx, delimiterIdx) : str.substring(tokenStartIdx) );
            result.add(token);
        } while( delimiterIdx != -1 );

        return result;
    }


    /**
     * join an array of strings given a seperator
     * @param separator the string to insert between each array element
     * @param strings the array of strings
     * @return a string, which is the joining of all array values with the separator
     */
    public static String join(String separator, String[] strings) {
        return join(separator, strings, 0, strings.length);
    }

    public static String join(String separator, String[] strings, int start, int end) {
        if ((end - start) == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[start]);
        for (int i = start + 1; i < end; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }

    /**
     * Returns a string of the form elt1.toString() [sep elt2.toString() ... sep elt.toString()] for a collection of
     * elti objects (note there's no actual space between sep and the elti elements).  Returns
     * "" if collection is empty.  If collection contains just elt, then returns elt.toString()
     *
     * @param separator the string to use to separate objects
     * @param objects a collection of objects.  the element order is defined by the iterator over objects
     * @param <T> the type of the objects
     * @return a non-null string
     */
    public static <T> String join(final String separator, final Collection<T> objects) {
        if (objects.isEmpty()) { // fast path for empty collection
            return "";
        } else {
            final Iterator<T> iter = objects.iterator();
            final T first = iter.next();

            if ( ! iter.hasNext() ) // fast path for singleton collections
                return first.toString();
            else { // full path for 2+ collection that actually need a join
                final StringBuilder ret = new StringBuilder(first.toString());
                while(iter.hasNext()) {
                    ret.append(separator);
                    ret.append(iter.next().toString());
                }
                return ret.toString();
            }
        }
    }

    public static String dupString(char c, int nCopies) {
        char[] chars = new char[nCopies];
        Arrays.fill(chars, c);
        return new String(chars);
    }

    public static byte[] dupBytes(byte b, int nCopies) {
        byte[] bytes = new byte[nCopies];
        Arrays.fill(bytes, b);
        return bytes;
    }

    // trim a string for the given character (i.e. not just whitespace)
    public static String trim(String str, char ch) {
        char[] array = str.toCharArray();


        int start = 0;
        while ( start < array.length && array[start] == ch )
            start++;

        int end = array.length - 1;
        while ( end > start && array[end] == ch )
            end--;

        return str.substring(start, end+1);
    }

    public static byte listMaxByte(List<Byte> quals) {
        if (quals.size() == 0) return 0;
        byte m = quals.get(0);
        for (byte b : quals) {
            m = b > m ? b : m;
        }
        return m;
    }


    // returns the maximum value in the array
    public static double findMaxEntry(double[] array) {
        return findIndexAndMaxEntry(array).first;
    }

    // returns the index of the maximum value in the array
    public static int findIndexOfMaxEntry(double[] array) {
        return findIndexAndMaxEntry(array).second;
    }

    // returns the the maximum value and its index in the array
    private static Pair<Double, Integer> findIndexAndMaxEntry(double[] array) {
        if ( array.length == 0 )
            return new Pair<Double, Integer>(0.0, -1);
        int index = 0;
        double max = array[0];
        for (int i = 1; i < array.length; i++) {
            if ( array[i] > max ) {
                max = array[i];
                index = i;
            }
        }
        return new Pair<Double, Integer>(max, index);
    }

    /**
     * Splits expressions in command args by spaces and returns the array of expressions.
     * Expressions may use single or double quotes to group any individual expression, but not both.
     * @param args Arguments to parse.
     * @return Parsed expressions.
     */
    public static String[] escapeExpressions(String args) {
        // special case for ' and " so we can allow expressions
        if (args.indexOf('\'') != -1)
            return escapeExpressions(args, "'");
        else if (args.indexOf('\"') != -1)
            return escapeExpressions(args, "\"");
        else
            return args.trim().split(" +");
    }

    /**
     * Splits expressions in command args by spaces and the supplied delimiter and returns the array of expressions.
     * @param args Arguments to parse.
     * @param delimiter Delimiter for grouping expressions.
     * @return Parsed expressions.
     */
    private static String[] escapeExpressions(String args, String delimiter) {
        String[] command = {};
        String[] split = args.split(delimiter);
        String arg;
        for (int i = 0; i < split.length - 1; i += 2) {
            arg = split[i].trim();
            if (arg.length() > 0) // if the unescaped arg has a size
                command = Utils.concatArrays(command, arg.split(" +"));
            command = Utils.concatArrays(command, new String[]{split[i + 1]});
        }
        arg = split[split.length - 1].trim();
        if (split.length % 2 == 1) // if the command ends with a delimiter
            if (arg.length() > 0) // if the last unescaped arg has a size
                command = Utils.concatArrays(command, arg.split(" +"));
        return command;
    }

    /**
     * Concatenates two String arrays.
     * @param A First array.
     * @param B Second array.
     * @return Concatenation of A then B.
     */
    public static String[] concatArrays(String[] A, String[] B) {
       String[] C = new String[A.length + B.length];
       System.arraycopy(A, 0, C, 0, A.length);
       System.arraycopy(B, 0, C, A.length, B.length);
       return C;
    }

    /**
     * Appends String(s) B to array A.
     * @param A First array.
     * @param B Strings to append.
     * @return A with B(s) appended.
     */
    public static String[] appendArray(String[] A, String... B) {
        return concatArrays(A, B);
    }

    /**
     * Returns indices of all occurrences of the specified symbol in the string
     * @param s Search string
     * @param ch Character to search for
     * @return Indices of all occurrences of the specified symbol
     */
    public static int[] indexOfAll(String s, int ch) {
        int[] pos = new int[64];
        int z = 0;

        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == ch) pos[z++] = i;
        }
        return reallocate(pos, z);
    }

    /**
     * Returns new (reallocated) integer array of the specified size, with content
     * of the original array <code>orig</code> copied into it. If <code>newSize</code> is
     * less than the size of the original array, only first <code>newSize</code> elements will be copied.
     * If new size is greater than the size of the original array, the content of the original array will be padded
     * with zeros up to the new size. Finally, if new size is the same as original size, no memory reallocation
     * will be performed and the original array will be returned instead.
     *
     * @param orig Original size.
     * @param newSize New Size.
     *
     * @return New array with length equal to newSize.
     */
    public static int[] reallocate(int[] orig, int newSize) {
        if (orig.length == newSize) return orig;
        int[] new_array = new int[newSize];
        int L = (newSize > orig.length ? orig.length : newSize);
        for (int i = 0; i < L; i++) new_array[i] = orig[i];
        return new_array;
    }


    /**
     * Returns a copy of array a, extended with additional n elements to the right (if n > 0 ) or -n elements to the
     * left (if n<0), copying the values form the original array. Newly added elements are filled with value v. Note that
     * if array a is being padded to the left, first (-n) elements of the returned array are v's, followed by the content of
     * array a.
     * @param a original array
     * @param n number of (v-filled) elements to append to a on the right (n>0) or on the left (n<0)
     * @param v element value
     * @return the extended copy of array a with additional n elements
     */
    public static byte [] extend(final byte[] a, int n, byte v) {

        byte [] newA;

        if ( n > 0 ) {
            newA = Arrays.copyOf(a, a.length+n);
            if ( v != 0) { // java pads with 0's for us, so there is nothing to do if v==0
                for ( int i = a.length; i < newA.length ; i++ ) newA[i] = v;
            }
            return newA;
        }

        // we are here only if n < 0:
        n = (-n);
        newA = new byte[ a.length + n ];
        int i;
        if ( v!= 0 ) {
            i = 0;
            for( ; i < n; i++ ) newA[i] = v;
        } else {
            i = n;
        }
        for ( int j = 0 ; j < a.length ; i++, j++) newA[i]=a[j];
        return newA;
    }


    /**
     * Returns a copy of array a, extended with additional n elements to the right (if n > 0 ) or -n elements to the
     * left (if n<0), copying the values form the original array. Newly added elements are filled with value v. Note that
     * if array a is padded to the left, first (-n) elements of the returned array are v's, followed by the content of
     * array a.
     * @param a original array
     * @param n number of (v-filled) elements to append to a on the right (n>0) or on the left (n<0)
     * @param v element value
     * @return the extended copy of array a with additional n elements
     */
    public static short [] extend(final short[] a, int n, short v) {

        short [] newA;

        if ( n > 0 ) {
            newA = Arrays.copyOf(a, a.length+n);
            if ( v != 0) { // java pads with 0's for us, so there is nothing to do if v==0
                for ( int i = a.length; i < newA.length ; i++ ) newA[i] = v;
            }
            return newA;
        }

        // we are here only if n < 0:
        n = (-n);
        newA = new short[ a.length + n ];
        int i;
        if ( v!= 0 ) {
            i = 0;
            for( ; i < n; i++ ) newA[i] = v;
        } else {
            i = n;
        }
        for ( int j = 0 ; j < a.length ; i++, j++) newA[i]=a[j];
        return newA;
    }

    /* TEST ME
        public static void main(String[] argv) {
            List<Integer> l1 = new LinkedList<Integer>();
            List<Integer> l2 = new ArrayList<Integer>();

            l1.add(1);
            l1.add(5);
            l1.add(3);
            l1.add(10);
            l1.add(4);
            l1.add(2);
            l2.add(1);
            l2.add(5);
            l2.add(3);
            l2.add(10);
            l2.add(4);
            l2.add(2);

            Predicate<Integer> p = new Predicate<Integer>() {
                public boolean apply(Integer i) {
                    return i > 2;
                }
            };
            filterInPlace(p, l1);
            filterInPlace(p, l2);

            for ( int i = 0 ; i < l1.size(); i++ ) System.out.print(" "+l1.get(i));
            System.out.println();
            for ( int i = 0 ; i < l2.size(); i++ ) System.out.print(" " + l2.get(i));
            System.out.println();

        }

    */

    /**
     * a helper method. Turns a single character string into a char.
     *
     * @param str the string
     *
     * @return a char
     */
    public static char stringToChar(String str) {
        if (str.length() != 1) throw new IllegalArgumentException("String length must be one");
        return str.charAt(0);
    }

    public static <T extends Comparable<T>> List<T> sorted(Collection<T> c) {
        return sorted(c, false);
    }

    public static <T extends Comparable<T>> List<T> sorted(Collection<T> c, boolean reverse) {
        List<T> l = new ArrayList<T>(c);
        Collections.sort(l);
        if ( reverse ) Collections.reverse(l);
        return l;
    }

    public static <T extends Comparable<T>, V> List<V> sorted(Map<T,V> c) {
        return sorted(c, false);
    }

    public static <T extends Comparable<T>, V> List<V> sorted(Map<T,V> c, boolean reverse) {
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);
        if ( reverse ) Collections.reverse(t);

        List<V> l = new ArrayList<V>();
        for ( T k : t ) {
            l.add(c.get(k));
        }
        return l;
    }

    public static <T extends Comparable<T>, V> String sortedString(Map<T,V> c) {
        List<T> t = new ArrayList<T>(c.keySet());
        Collections.sort(t);

        List<V> l = new ArrayList<V>();
        List<String> pairs = new ArrayList<String>();
        for ( T k : t ) {
            pairs.add(k + "=" + c.get(k));
        }

        return "{" + join(", ", pairs) + "}";
    }

    /**
     * Reverse a byte array of bases
     *
     * @param bases  the byte array of bases
     * @return the reverse of the base byte array
     */
    static public byte[] reverse(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    static public final <T> List<T> reverse(final List<T> l) {
        final List<T> newL = new ArrayList<T>(l);
        Collections.reverse(newL);
        return newL;
    }

    /**
     * Reverse an int array of bases
     *
     * @param bases  the int array of bases
     * @return the reverse of the base int array
     */
    static public int[] reverse(int[] bases) {
        int[] rcbases = new int[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    /**
     * Reverse (NOT reverse-complement!!) a string
     *
     * @param bases  input string
     * @return the reversed string
     */
    static public String reverse(String bases) {
        return new String( reverse( bases.getBytes() )) ;
    }

    public static byte[] charSeq2byteSeq(char[] seqIn) {
        byte[] seqOut = new byte[seqIn.length];
        for ( int i = 0; i < seqIn.length; i++ ) {
            seqOut[i] = (byte)seqIn[i];
        }
        return seqOut;
    }

    public static boolean isFlagSet(int value, int flag) {
        return ((value & flag) == flag);
    }

    /**
     * Helper utility that calls into the InetAddress system to resolve the hostname.  If this fails,
     * unresolvable gets returned instead.
     *
     * @return
     */
    public static final String resolveHostname() {
        try {
            return InetAddress.getLocalHost().getCanonicalHostName();
        }
        catch (java.net.UnknownHostException uhe) { // [beware typo in code sample -dmw]
            return "unresolvable";
            // handle exception
        }
    }
}
