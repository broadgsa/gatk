package org.broadinstitute.sting.utils;

import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;
import org.apache.log4j.Logger;

import java.io.File;
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

    public static void warnUser(final String msg) {
        logger.warn(String.format("********************************************************************************"));
        logger.warn(String.format("* WARNING:"));
        logger.warn(String.format("*"));
        prettyPrintWarningMessage(msg);
        logger.warn(String.format("********************************************************************************"));
    }

    public static void scareUser(final String msg) {
        //System.out.printf("********************************************************************************%n");
        //System.out.printf("* ERROR:%n");
        //System.out.printf("*%n");
        //System.out.printf("* %s%n", msg);
        //System.out.printf("********************************************************************************%n");
        logger.fatal(msg);
        throw new StingException(msg);
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

    /**
     * pretty print the warning message supplied
     *
     * @param message the message
     */
    private static void prettyPrintWarningMessage(String message) {
        StringBuilder builder = new StringBuilder(message);
        while (builder.length() > 70) {
            int space = builder.lastIndexOf(" ", 70);
            if (space <= 0) space = 70;
            logger.warn(String.format("* %s", builder.substring(0, space)));
            builder.delete(0, space + 1);
        }
        logger.warn(String.format("* %s", builder));
    }

    public static SAMFileHeader copySAMFileHeader(SAMFileHeader toCopy) {
        SAMFileHeader copy = new SAMFileHeader();

        copy.setSortOrder(toCopy.getSortOrder());
        copy.setGroupOrder(toCopy.getGroupOrder());
        copy.setProgramRecords(toCopy.getProgramRecords());
        copy.setReadGroups(toCopy.getReadGroups());
        copy.setSequenceDictionary(toCopy.getSequenceDictionary());

        for (Map.Entry<String, Object> e : toCopy.getAttributes())
            copy.setAttribute(e.getKey(), e.getValue());

        return copy;
    }

    public static SAMFileWriter createSAMFileWriterWithCompression(SAMFileHeader header, boolean presorted, String file, int compression) {
        if (file.endsWith(".bam"))
            return new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(file), compression);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, new File(file));
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

    public static boolean isPlatformRead(SAMRecord read, String name) {
        SAMReadGroupRecord readGroup = read.getReadGroup();
        if (readGroup != null) {
            Object readPlatformAttr = readGroup.getAttribute("PL");
            if (readPlatformAttr != null)
                return readPlatformAttr.toString().toUpperCase().contains(name);
        }
        return false;
    }


    public static boolean is454Read(SAMRecord read) {
        return isPlatformRead(read, "454");
    }

    public static boolean isSOLiDRead(SAMRecord read) {
        return isPlatformRead(read, "SOLID");
    }

    public static boolean isSLXRead(SAMRecord read) {
        return isPlatformRead(read, "ILLUMINA");
    }

    private static final Map<Integer, String> readFlagNames
            = new HashMap<Integer, String>();

    static {
        readFlagNames.put(0x1, "Paired");
        readFlagNames.put(0x2, "Proper");
        readFlagNames.put(0x4, "Unmapped");
        readFlagNames.put(0x8, "MateUnmapped");
        readFlagNames.put(0x10, "Forward");
        //readFlagNames.put(0x20, "MateForward");
        readFlagNames.put(0x4, "FirstOfPair");
        readFlagNames.put(0x8, "SecondOfPair");
        readFlagNames.put(0x100, "NotPrimary");
        readFlagNames.put(0x200, "NON-PF");
        readFlagNames.put(0x400, "Duplicate");
    }

    public static String readFlagsAsString(SAMRecord rec) {
        String flags = "";
        for (int flag : readFlagNames.keySet()) {
            if ((rec.getFlags() & flag) != 0) {
                flags += readFlagNames.get(flag) + " ";
            }
        }
        return flags;
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

    //public static String join(String separator, Collection<String> strings) {
    //    return join( separator, strings.toArray(new String[0]) );
    //}

    public static <T> String join(String separator, Collection<T> objects) {
        ArrayList<String> strs = new ArrayList<String>();
        for (Object x : objects)
            strs.add(x.toString());
        return join(separator, strs.toArray(new String[0]));
    }

    public static double average(List<Long> vals, int maxI) {
        long sum = 0L;

        int i = 0;
        for (long x : vals) {
            if (i > maxI)
                break;
            sum += x;
            i++;
            //System.out.printf(" %d/%d", sum, i);
        }

        //System.out.printf("Sum = %d, n = %d, maxI = %d, avg = %f%n", sum, i, maxI, (1.0 * sum) / i);

        return (1.0 * sum) / i;
    }

    public static double averageDouble(List<Double> vals, int maxI) {
        double sum = 0.0;

        int i = 0;
        for (double x : vals) {
            if (i > maxI)
                break;
            sum += x;
            i++;
        }
        return (1.0 * sum) / i;
    }

    public static double average(List<Long> vals) {
        return average(vals, vals.size());
    }

    public static double averageDouble(List<Double> vals) {
        return averageDouble(vals, vals.size());
    }

    // Java Generics can't do primitive types, so I had to do this the simplistic way

    public static Integer[] SortPermutation(final int[] A) {
        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                if (A[a.intValue()] < A[b.intValue()]) {
                    return -1;
                }
                if (A[a.intValue()] == A[b.intValue()]) {
                    return 0;
                }
                if (A[a.intValue()] > A[b.intValue()]) {
                    return 1;
                }
                return 0;
            }
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static Integer[] SortPermutation(final double[] A) {
        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                if (A[a.intValue()] < A[b.intValue()]) {
                    return -1;
                }
                if (A[a.intValue()] == A[b.intValue()]) {
                    return 0;
                }
                if (A[a.intValue()] > A[b.intValue()]) {
                    return 1;
                }
                return 0;
            }
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static <T extends Comparable> Integer[] SortPermutation(List<T> A) {
        final Object[] data = A.toArray();

        class comparator implements Comparator<Integer> {
            public int compare(Integer a, Integer b) {
                return ((T) data[a]).compareTo(data[b]);
            }
        }
        Integer[] permutation = new Integer[A.size()];
        for (int i = 0; i < A.size(); i++) {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }


    public static int[] PermuteArray(int[] array, Integer[] permutation) {
        int[] output = new int[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static double[] PermuteArray(double[] array, Integer[] permutation) {
        double[] output = new double[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static Object[] PermuteArray(Object[] array, Integer[] permutation) {
        Object[] output = new Object[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static String[] PermuteArray(String[] array, Integer[] permutation) {
        String[] output = new String[array.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static <T> List<T> PermuteList(List<T> list, Integer[] permutation) {
        List<T> output = new ArrayList<T>();
        for (int i = 0; i < permutation.length; i++) {
            output.add(list.get(permutation[i]));
        }
        return output;
    }


    /** Draw N random elements from list. */
    public static <T> List<T> RandomSubset(List<T> list, int N) {
        if (list.size() <= N) {
            return list;
        }

        java.util.Random random = new java.util.Random();

        int idx[] = new int[list.size()];
        for (int i = 0; i < list.size(); i++) {
            idx[i] = random.nextInt();
        }

        Integer[] perm = SortPermutation(idx);

        List<T> ans = new ArrayList<T>();
        for (int i = 0; i < N; i++) {
            ans.add(list.get(perm[i]));
        }

        return ans;
    }

    // lifted from the internet 
    // http://www.cs.princeton.edu/introcs/91float/Gamma.java.html
    public static double logGamma(double x) {
        double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
        double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                + 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
        return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }

    public static double percentage(double x, double base) {
        return (base > 0 ? (x / base) * 100.0 : 0);
    }

    public static double percentage(int x, int base) {
        return (base > 0 ? ((double) x / (double) base) * 100.0 : 0);
    }

    public static double percentage(long x, long base) {
        return (base > 0 ? ((double) x / (double) base) * 100.0 : 0);
    }

    public static String dupString(char c, int nCopies) {
        char[] chars = new char[nCopies];
        Arrays.fill(chars, c);
        return new String(chars);
    }

    public static int countOccurrences(char c, String s) {
        int count = 0;
        for (int i = 0; i < s.length(); i++) {
            count += s.charAt(i) == c ? 1 : 0;
        }
        return count;
    }

    public static <T> int countOccurrences(T x, List<T> l) {
        int count = 0;
        for (T y : l) {
            if (x.equals(y)) count++;
        }

        return count;
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

    public static String[] concatArrays(String[] A, String[] B) {
       String[] C = new String[A.length + B.length];
       System.arraycopy(A, 0, C, 0, A.length);
       System.arraycopy(B, 0, C, A.length, B.length);
       return C;
    }


    /** Returns indices of all occurrences of the specified symbol in the string */
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
     * @param orig
     * @param newSize
     *
     * @return
     */
    public static int[] reallocate(int[] orig, int newSize) {
        if (orig.length == newSize) return orig;
        int[] new_array = new int[newSize];
        int L = (newSize > orig.length ? orig.length : newSize);
        for (int i = 0; i < L; i++) new_array[i] = orig[i];
        return new_array;
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
}




