package org.broadinstitute.sting.utils;

import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequenceFile;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;

import org.apache.log4j.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:12:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class Utils {
    /**
     * our log, which we want to capture anything from this class
     */
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
        throw new RuntimeException(msg);
    }

    public static <T> List<T> cons(final T elt, final List<T> l) {
        List<T> l2 = new ArrayList<T>();
        l2.add(elt);
        if ( l != null ) l2.addAll(l);
        return l2;
    }

    /**
     * pretty print the warning message supplied
     * @param message the message
     */
    private static void prettyPrintWarningMessage(String message) {
        StringBuilder builder = new StringBuilder(message);
        while (builder.length() > 70) {
            int space = builder.lastIndexOf(" ", 70);
            if (space <= 0) space = 70;
            logger.warn(String.format("* %s", builder.substring(0,space)));
            builder.delete(0,space + 1);
        }
        logger.warn(String.format("* %s", builder));
    }

    public static SAMFileHeader copySAMFileHeader( SAMFileHeader toCopy ) {
        SAMFileHeader copy = new SAMFileHeader();

        copy.setSortOrder(toCopy.getSortOrder());
        copy.setGroupOrder(toCopy.getGroupOrder());
        copy.setProgramRecords(toCopy.getProgramRecords());
        copy.setReadGroups(toCopy.getReadGroups());
        copy.setSequenceDictionary(toCopy.getSequenceDictionary());

        for ( Map.Entry<String, Object> e : toCopy.getAttributes())
            copy.setAttribute(e.getKey(), e.getValue());
        
        return copy;
    }

    public static SAMFileWriter createSAMFileWriterWithCompression(SAMFileHeader header, boolean presorted, String file, int compression) {
        if (file.endsWith(".bam"))
            return new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(file), compression);
        return new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, new File(file));
    }

    /**
     * Returns a new list built from those objects found in collection <c> that satisfy the
     * predicate ( i.e. pred.apply() is true for the objects in th eresulting list ).
     *
     * @param pred filtering condition ( objects, for which pred.apply() is true pass the filter )
     * @param c    collection to filter (will not be modified)
     * @return new list built from elements of <c> passing the filter
     * @see #filterInPlace(Predicate pred, Collection c)
     */
    public static <T> List<T> filter(Predicate<? super T> pred, Collection<T> c) {
        List<T> filtered = new ArrayList<T>();
        // loop through all the elements in c
        for (T obj : c) {
            // if the predicate is true for the current element
            if (pred.apply(obj)) {
                // append it to the result list
                filtered.add(obj);
            }
        }
        return filtered;
    }

    /**
     * Removes from the collection <c> all the elements that do not pass the filter (i.e. those elements,
     * for which pred.apply() is false ). This is an in-place method - the argument is modified, and no new
     * objects are created/copied. Collection's iterator (as returned by iterator()) must implement
     * optional remove() interface method that allows multiple subsequent removals of elements from the
     * underlying collection (this is the standard contract). This method
     * works best for collections that support cheap, constant time
     * object removal (such as LinkedList, HashSet etc.). It is also specifically designed to
     * detect ArrayLists and use optimized strategy for them. However
     * with other, custom lists that 1) do not inherit (are not instanceof) from ArrayList and 2) do not implement
     * fast (constant time) remove() operation, the performance can degrade significantly (linear traversal times,
     * e.g., linear removal ~ N^2).
     *
     * @param pred filtering condition (only elements, for which pred.apply() is true will be kept in the collection)
     * @param c    collection to filter (will be modified - should be mutable and should implement remove() )
     * @return reference to the same (modified) collection <c>
     * @see #filter(Predicate pred, Collection c)
     */
    public static <T> Collection<T> filterInPlace(Predicate<? super T> pred, Collection<T> c) {
        if (c instanceof ArrayList) {
            // arraylists are a special case that we know how to process efficiently
            // (generic implementation below removes one element at a time and is not well suited
            // for ArrayLists
            List<T> list = (List<T>) c;
            int j = 0; // copy-to location
            // perform one linear pass copying forward all elements that pass the filter,
            // so that the head of the list is continuous sequence of such elements:
            for (int i = 0; i < list.size(); i++) {
                // if object passes, copy it forward and increment j (=copy-to location);
                // otherwise keep the same copy-to location and move on to the next element
                if (pred.apply(list.get(i))) list.set(j++, list.get(i));
            }
            // j now points to first unused copy-to location; elements 0...j-1 pass the filter
            list.subList(j, list.size()).clear(); // remove tail of the list
        }
/*
        // loop through all the elements in c
        for (T obj : c) {
            // if the predicate is false for the current element
            if (! pred.apply(obj)) {
                // remove that element from the collection
                c.remove(obj);
            }
        }
 */
        Iterator<T> it = c.iterator();
        while (it.hasNext()) {
            if (pred.apply(it.next())) continue;
            it.remove();
        }
        return c;
    }

    public static ArrayList<Byte> subseq(char[] fullArray) {
        byte[] fullByteArray = new byte[fullArray.length];
        StringUtil.charsToBytes(fullArray,0,fullArray.length,fullByteArray,0);
        return subseq(fullByteArray);
    }

    public static ArrayList<Byte> subseq(byte[] fullArray) {
        return subseq(fullArray, 0, fullArray.length-1);
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

    public static boolean is454Read(SAMRecord read, SAMFileHeader header) {
        Object readGroupAttr = read.getAttribute("RG");
        if ( readGroupAttr != null ) {
            SAMReadGroupRecord readGroup = header.getReadGroup(readGroupAttr.toString());
            if ( readGroup != null ) {
            Object readPlatformAttr = readGroup.getAttribute("PL");
            if ( readPlatformAttr != null )
                return readPlatformAttr.toString().toUpperCase().contains("454");
            }
        }
        return false;
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

    public static String join(String separator, String[] strings) {
        return join(separator, strings, 0, strings.length);
    }

    public static String join(String separator, String[] strings, int start, int end) {
        if ((end - start) == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[start]);
        for (int i = start+1; i < end; ++i) {
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
                if (A[a.intValue()] < A[ b.intValue() ]) {
                    return -1;
                }
                if (A[ a.intValue() ] == A[ b.intValue() ]) {
                    return 0;
                }
                if (A[ a.intValue() ] > A[ b.intValue() ]) {
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

    public static <T> List<T> PermuteList(List<T> list, Integer[] permutation) 
    {
        List<T> output = new ArrayList<T>();
        for (int i = 0; i < permutation.length; i++) {
            output.add(list.get(permutation[i]));
        }
        return output;
    }


    /** Draw N random elements from list. */
    public static <T> List<T> RandomSubset(List<T> list, int N)
    {
        if (list.size() <= N) { return list; }

        java.util.Random random = new java.util.Random();

        int idx[] = new int[list.size()];
        for (int i = 0; i < list.size(); i++) { idx[i] = random.nextInt(); }

        Integer[] perm = SortPermutation(idx);        

        List<T> ans = new ArrayList<T>();
        for (int i = 0; i < N; i++) { ans.add(list.get(perm[i])); }

        return ans;
    }

    // lifted from the internet 
    // http://www.cs.princeton.edu/introcs/91float/Gamma.java.html
    public static double logGamma(double x) 
    {
		double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
		double ser = 1.0 + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
		                 + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
		                 +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5);
		return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
    }

    public static double percentage(double x, double base) { return (base> 0 ? (x/base)*100.0 : 0); }
    public static double percentage(int x, int base) { return (base> 0 ? ((double)x/(double)base)*100.0 : 0); }
    public static double percentage(long x, long base) { return (base> 0 ? ((double)x/(double)base)*100.0 : 0); }

    public static String dupString( char c, int nCopies ) {
        char[] chars = new char[nCopies];
        for ( int i = 0; i < nCopies; i++ ) chars[i] = c;
        //System.out.printf("chars is %s%n", new String(chars));
        return new String(chars);
    }

    public static int countOccurances(char c, String s) {
        int count = 0;
        for ( int i = 0; i < s.length(); i++ ) {
            count += s.charAt(i) == c ? 1 : 0;
        }
        return count;
    }

    public static byte listMaxByte(List<Byte> quals) {
        if ( quals.size() == 0 ) return 0;
        byte m = quals.get(0);
        for ( byte b : quals ) {
            m = b > m ? b : m;
        }
        return m;
    }


    /** Returns indices of all occurrences of the specified symbol in the string */
    public static int[] indexOfAll(String s, int ch) {
    	int[] pos = new int[64];
    	int z = 0;
    	
    	for ( int i = 0 ; i < s.length() ; i++ ) {
    		if ( s.charAt(i) == ch ) pos[z++] = i; 
    	}
    	return reallocate(pos,z);
    }
    
    /** Returns new (reallocated) integer array of the specified size, with content
     * of the original array <code>orig</code> copied into it. If <code>newSize</code> is
     * less than the size of the original array, only first <code>newSize</code> elements will be copied.
     * If new size is greater than the size of the original array, the content of the original array will be padded
     * with zeros up to the new size. Finally, if new size is the same as original size, no memory reallocation 
     * will be performed and the original array will be returned instead. 
     * @param orig
     * @param newSize
     * @return
     */
    public static int[] reallocate(int[] orig, int newSize) {
    	if ( orig.length == newSize ) return orig;
    	int[] new_array = new int[newSize];
    	int L = ( newSize > orig.length ? orig.length : newSize );
    	for ( int i = 0 ; i < L ; i++ ) new_array[i] = orig[i];
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
}




