package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.io.File;

/**
 * Parse text representations of interval strings that
 * can appear in Sting-based applications.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalUtils {
    /**
     * Turns a set of strings describing intervals into a parsed set of intervals.  Valid string elements can be files,
     * intervals in samtools notation (chrA:B-C), or some combination of the above separated by semicolons.  Additionally,
     * 'all' can be supplied to indicate all possible intervals, but 'all' must be exclusive of all other interval
     * specifications.
     *
     * @param argList A list of strings containing interval data.
     * @return an unsorted, unmerged representation of the given intervals.  Null is used to indicate that all intervals should be used. 
     */
    public static GenomeLocSortedSet parseIntervalArguments(List<String> argList) {
        List<GenomeLoc> rawIntervals = new ArrayList<GenomeLoc>();    // running list of raw GenomeLocs

        if (argList != null) { // now that we can be in this function if only the ROD-to-Intervals was provided, we need to
                               // ensure that the arg list isn't null before looping.
            for (String argument : argList) {

                // if any interval argument is '-L all', consider all loci by returning no intervals
                if (argument.equals("all")) {
                    if (argList.size() != 1) {
                        // throw error if '-L all' is not only interval - potentially conflicting commands
                        throw new StingException(String.format("Conflicting arguments: Intervals given along with \"-L all\""));
                    }
                    return null;
                }

                // separate argument on semicolon first
                for (String fileOrInterval : argument.split(";")) {

                    // if it's a file, add items to raw interval list
                    if (isFile(fileOrInterval))
                        rawIntervals.addAll(GenomeLocParser.intervalFileToList(fileOrInterval));

                        // otherwise treat as an interval -> parse and add to raw interval list
                    else {
                        rawIntervals.add(GenomeLocParser.parseGenomeInterval(fileOrInterval));
                    }
                }
            }
        }

        return GenomeLocSortedSet.createSetFromList(rawIntervals);
    }

    /**
     * Sorts and merges an interval list.  Multiple techniques are available for merging: ALL, which combines
     * all overlapping and abutting intervals into an interval that spans the union of all covered bases, and
     * OVERLAPPING_ONLY, which unions overlapping intervals but keeps abutting intervals separate.
     *
     * @param intervals A collection of intervals to merge.
     * @param mergingRule A descriptor for the type of merging to perform.
     * @return A sorted, merged version of the intervals passed in.
     */
    public static GenomeLocSortedSet sortAndMergeIntervals(GenomeLocSortedSet intervals, IntervalMergingRule mergingRule) {
        List<GenomeLoc> intervalList = intervals.toList();

        // sort raw interval list
        Collections.sort(intervalList);
        // now merge raw interval list
        intervalList = GenomeLocParser.mergeIntervalLocations(intervalList, mergingRule);

        return GenomeLocSortedSet.createSetFromList(intervalList);
    }

    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @return true if the token looks like a filename, or false otherwise.
     */
    private static boolean isFile(String str) {
        // should we define list of file extensions as a public array somewhere?
        // is regex or endsiwth better?
        if (str.toUpperCase().endsWith(".BED") || str.toUpperCase().endsWith(".LIST") ||
                str.toUpperCase().endsWith(".PICARD") || str.toUpperCase().endsWith(".INTERVAL_LIST")
                || str.toUpperCase().endsWith(".INTERVALS"))
            return true;

        if(new File(str).exists())
            throw new StingException("Interval argument looks like a filename, but does not have one of " +
                                     "the supported extensions (.bed, .list, .picard, .interval_list, or .intervals).  " +
                                     "Please rename your file with the appropriate extension.");

        else return false;
    }    
}
