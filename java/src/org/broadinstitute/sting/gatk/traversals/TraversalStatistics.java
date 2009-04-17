package org.broadinstitute.sting.gatk.traversals;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 8, 2009
 * Time: 4:13:40 PM
 *
 * Holds a bunch of basic information about the traversal.
 * TODO: Make this a class that can be passed around from the TraversalEngine to other entries that want to update it.
 */
public class TraversalStatistics {
    // Number of records (loci, reads) we've processed
    public static long nRecords;
    // How many reads have we processed, along with those skipped for various reasons
    public static int nReads;
    public static int nSkippedReads;
    public static int nUnmappedReads;
    public static int nNotPrimary;
    public static int nBadAlignments;
    public static int nSkippedIndels;
    public static int nDuplicates;

    static {
        reset();
    }

    public static void reset() {
        nRecords = 0;
        nReads = 0;
        nSkippedReads = 0;
        nUnmappedReads = 0;
        nNotPrimary = 0;
        nBadAlignments = 0;
        nSkippedIndels = 0;
        nDuplicates = 0;
    }
}
