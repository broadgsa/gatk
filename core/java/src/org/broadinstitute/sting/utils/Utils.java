package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Feb 24, 2009
 * Time: 10:12:31 AM
 * To change this template use File | Settings | File Templates.
 */
public class Utils {
    public static void warnUser(final String msg) {
        System.out.printf("********************************************************************************%n");
        System.out.printf("* WARNING:%n");
        System.out.printf("*%n");
        System.out.printf("* %s%n", msg);
        System.out.printf("********************************************************************************%n");                
    }

    public static void scareUser(final String msg) {
        System.out.printf("********************************************************************************%n");
        System.out.printf("* ERROR:%n");
        System.out.printf("*%n");
        System.out.printf("* %s%n", msg);
        System.out.printf("********************************************************************************%n");
        throw new RuntimeException(msg);
    }

    /** Returns a new list built from those objects found in collection <c> that satisfy the
     * predicate ( i.e. pred.apply() is true for the objects in th eresulting list ).
      * @param pred filtering condition ( objects, for which pred.apply() is true pass the filter )
     * @param c collection to filter (will not be modified)
     * @return  new list built from elements of <c> passing the filter
     * @see #filterInPlace(Predicate pred, Collection c)
     */
    public static <T> List<T> filter(Predicate pred, Collection<T> c) {
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

    /** Removes from the collection <c> all the elements that do not pass the filter (i.e. those elements,
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
     * @param pred filtering condition (only elements, for which pred.apply() is true will be kept in the collection)
     * @param c collection to filter (will be modified - should be mutable and should implement remove() )
     * @return reference to the same (modified) collection <c>
     * @see #filter(Predicate pred, Collection c)
     */
    public static <T> Collection<T> filterInPlace(Predicate pred, Collection<T> c) {
        if ( c instanceof ArrayList ) {
            // arraylists are a special case that we know how to process efficiently
            // (generic implementation below removes one element at a time and is not well suited
            // for ArrayLists
            List<T> list = (List<T>)c;
            int j = 0; // copy-to location
            // perform one linear pass copying forward all elements that pass the filter,
            // so that the head of the list is continuous sequence of such elements:
            for ( int i = 0 ; i < list.size() ; i++ ) {
                // if object passes, copy it forward and increment j (=copy-to location);
                // otherwise keep the same copy-to location and move on to the next element
                 if ( pred.apply(list.get(i)) ) list.set(j++,list.get(i));
            }
            // j now points to first unused copy-to location; elements 0...j-1 pass the filter
            list.subList(j,list.size()).clear(); // remove tail of the list
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
        while ( it.hasNext() ) {
            if ( pred.apply(it.next() )) continue;
            it.remove();
        }
        return c;
    }

    public static ArrayList<Byte> subseq(byte[] fullArray) {
        return subseq(fullArray, 0, fullArray.length);
    }
    
    public static ArrayList<Byte> subseq(byte[] fullArray, int start, int end) {
        ArrayList<Byte> dest = new ArrayList<Byte>(end-start+1);
        for ( int i = start; i < end; i++ ) {
            dest.add(fullArray[i]);
        }
        return dest;
    }

    public static String baseList2string(List<Byte> bases) {
        byte[] basesAsbytes = new byte[bases.size()];
        int i = 0;
        for ( Byte b : bases ) {
            basesAsbytes[i] = b;
            i++;
        }
        return new String(basesAsbytes);
    }

    public static GenomeLoc genomicLocationOf( final SAMRecord read ) {
        return new GenomeLoc( read.getReferenceName(), read.getAlignmentStart() );
    }

    private static final Map<Integer,String> readFlagNames
          = new HashMap<Integer,String>();

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
        for ( int flag : readFlagNames.keySet() ) {
            if ( ( rec.getFlags() & flag ) != 0 ) {
                flags += readFlagNames.get(flag) + " ";
            }
        }
        return flags;
    }

    public static String join(String separator, String[] strings) {
        if (strings.length == 0) {
            return "";
        }
        StringBuilder ret = new StringBuilder(strings[0]);
        for (int i = 1; i < strings.length; ++i) {
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
        for ( Object x : objects )
            strs.add(x.toString());
        return join( separator, strs.toArray(new String[0]) );
    }

    public static double average(List<Long> vals, int maxI) {
        long sum = 0L;

        int i = 0;
        for ( long x : vals ) {
            if ( i > maxI )
                break;
            sum += x;
            i++;
            //System.out.printf(" %d/%d", sum, i);
        }

        //System.out.printf("Sum = %d, n = %d, maxI = %d, avg = %f%n", sum, i, maxI, (1.0 * sum) / i);

        return (1.0 * sum) / i;
    }

    public static double average(List<Long> vals) {
        return average(vals, vals.size());
    }

    public static boolean setupRefContigOrdering(final ReferenceSequenceFile refFile) {
        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary();
        HashMap<String, Integer> refContigOrdering = new HashMap<String, Integer>();

        if ( refContigs != null ) {
            int i = 0;
            System.out.printf("Prepared reference sequence contig dictionary%n  order ->");
            for ( SAMSequenceRecord contig : refContigs ) {
                System.out.printf(" %s", contig.getSequenceName());
                refContigOrdering.put(contig.getSequenceName(), i);
                i++;
            }
            System.out.printf("%n  Total elements -> %d%n", refContigOrdering.size());
        }
        
        GenomeLoc.setContigOrdering(refContigOrdering);
        return refContigs != null;
    }

    // Java Generics can't do primitive types, so I had to do this the simplistic way
    
    public static Integer[] SortPermutation(final int[] A)
    {
        class comparator implements Comparator
        { 
            public int compare(Object a, Object b)
            {
                if (A[(Integer)a] <  A[(Integer)b]) { return -1; }
                if (A[(Integer)a] == A[(Integer)b]) { return  0; }
                if (A[(Integer)a] >  A[(Integer)b]) { return  1; }
                return 0;
            } 
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++)
        {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static Integer[] SortPermutation(final double[] A)
    {
        class comparator implements Comparator
        { 
            public int compare(Object a, Object b)
            {
                if (A[(Integer)a] <  A[(Integer)b]) { return -1; }
                if (A[(Integer)a] == A[(Integer)b]) { return  0; }
                if (A[(Integer)a] >  A[(Integer)b]) { return  1; }
                return 0;
            } 
        }
        Integer[] permutation = new Integer[A.length];
        for (int i = 0; i < A.length; i++)
        {
            permutation[i] = i;
        }
        Arrays.sort(permutation, new comparator());
        return permutation;
    }

    public static int[] PermuteArray(int[] array, Integer[] permutation)
    {
        int[] output = new int[array.length];
        for (int i = 0; i < output.length; i++)
        {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static double[] PermuteArray(double[] array, Integer[] permutation)
    {
        double[] output = new double[array.length];
        for (int i = 0; i < output.length; i++)
        {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static Object[] PermuteArray(Object[] array, Integer[] permutation)
    {
        Object[] output = new Object[array.length];
        for (int i = 0; i < output.length; i++)
        {
            output[i] = array[permutation[i]];
        }
        return output;
    }

    public static String[] PermuteArray(String[] array, Integer[] permutation)
    {
        String[] output = new String[array.length];
        for (int i = 0; i < output.length; i++)
        {
            output[i] = array[permutation[i]];
        }
        return output;
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




