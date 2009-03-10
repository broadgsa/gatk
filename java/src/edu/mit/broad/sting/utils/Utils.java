package edu.mit.broad.sting.utils;

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

    public static String join(String separator, Collection<String> strings) {
        return join( separator, strings.toArray(new String[0]) );
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

    public static void setupRefContigOrdering(final ReferenceSequenceFile refFile) { 
        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary();
        HashMap<String, Integer> refContigOrdering = new HashMap<String, Integer>();

        int i = 0;
        System.out.printf("Prepared reference sequence contig dictionary%n  order ->");
        for ( SAMSequenceRecord contig : refContigs ) {
            System.out.printf(" %s", contig.getSequenceName());
            refContigOrdering.put(contig.getSequenceName(), i);
            i++;
        }
        System.out.printf("%n  Total elements -> %d%n", refContigOrdering.size());
        
        GenomeLoc.setContigOrdering(refContigOrdering);
    }
}
