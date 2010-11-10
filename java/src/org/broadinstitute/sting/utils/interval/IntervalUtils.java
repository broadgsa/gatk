package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.LinkedList;
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
     * @param allowEmptyIntervalList If false instead of an empty interval list will return null.
     * @return an unsorted, unmerged representation of the given intervals.  Null is used to indicate that all intervals should be used. 
     */
    public static List<GenomeLoc> parseIntervalArguments(GenomeLocParser parser, List<String> argList, boolean allowEmptyIntervalList) {
        List<GenomeLoc> rawIntervals = new ArrayList<GenomeLoc>();    // running list of raw GenomeLocs

        if (argList != null) { // now that we can be in this function if only the ROD-to-Intervals was provided, we need to
                               // ensure that the arg list isn't null before looping.
            for (String argument : argList) {

                // if any interval argument is '-L all', consider all loci by returning no intervals
                if (argument.equals("all")) {
                    if (argList.size() != 1) {
                        // throw error if '-L all' is not only interval - potentially conflicting commands
                        throw new UserException.CommandLineException(String.format("Conflicting arguments: Intervals given along with \"-L all\""));
                    }
                    return null;
                }

                // separate argument on semicolon first
                for (String fileOrInterval : argument.split(";")) {

                    // if it's a file, add items to raw interval list
                    if (isIntervalFile(fileOrInterval)) {
                        try {
                            rawIntervals.addAll(parser.intervalFileToList(fileOrInterval, allowEmptyIntervalList));
                        }
                        catch (Exception e) {
                            throw new UserException.MalformedFile(fileOrInterval, "Interval file could not be parsed in either format.", e);
                        }
                    }

                        // otherwise treat as an interval -> parse and add to raw interval list
                    else {
                        rawIntervals.add(parser.parseGenomeInterval(fileOrInterval));
                    }
                }
            }
        }

        return rawIntervals;
    }

    /**
     * merge two interval lists, using an interval set rule
     * @param setOne a list of genomeLocs, in order (cannot be NULL)
     * @param setTwo a list of genomeLocs, also in order (cannot be NULL)
     * @param rule the rule to use for merging, i.e. union, intersection, etc
     * @return a list, correctly merged using the specified rule
     */
    public static List<GenomeLoc> mergeListsBySetOperator(List<GenomeLoc> setOne, List<GenomeLoc> setTwo, IntervalSetRule rule) {
        // shortcut, if either set is zero, return the other set
        if (setOne == null || setOne.size() == 0 || setTwo == null || setTwo.size() == 0) return (setOne == null || setOne.size() == 0) ? setTwo : setOne;

        // if we're set to UNION, just add them all
        if (rule == IntervalSetRule.UNION) {
            setOne.addAll(setTwo);
            return setOne;
        }

        // else we're INTERSECTION, create two indexes into the lists
        int iOne = 0;
        int iTwo = 0;

        // our master list, since we can't guarantee removal time in a generic list
        LinkedList<GenomeLoc> retList = new LinkedList<GenomeLoc>();

        // merge the second into the first using the rule
        while (iTwo < setTwo.size() && iOne < setOne.size())
            // if the first list is ahead, drop items off the second until we overlap
            if (setTwo.get(iTwo).isBefore(setOne.get(iOne)))
                iTwo++;
            // if the second is ahead, drop intervals off the first until we overlap
            else if (setOne.get(iOne).isBefore(setTwo.get(iTwo)))
                iOne++;
            // we overlap, intersect the two intervals and add the result.  Then remove the interval that ends first.
            else {
                retList.add(setOne.get(iOne).intersect(setTwo.get(iTwo)));
                if (setOne.get(iOne).getStop() < setTwo.get(iTwo).getStop()) iOne++;
                else iTwo++;
            }
        
        // we don't need to add the rest of remaining locations, since we know they don't overlap. return what we have
        return retList;
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
    public static GenomeLocSortedSet sortAndMergeIntervals(GenomeLocParser parser, List<GenomeLoc> intervals, IntervalMergingRule mergingRule) {
        // sort raw interval list
        Collections.sort(intervals);
        // now merge raw interval list
        intervals = parser.mergeIntervalLocations(intervals, mergingRule);

        return GenomeLocSortedSet.createSetFromList(parser,intervals);
    }

    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(String str) {
        // should we define list of file extensions as a public array somewhere?
        // is regex or endsiwth better?
        File file = new File(str);
        if (str.toUpperCase().endsWith(".BED") || str.toUpperCase().endsWith(".LIST") ||
                str.toUpperCase().endsWith(".PICARD") || str.toUpperCase().endsWith(".INTERVAL_LIST")
                || str.toUpperCase().endsWith(".INTERVALS")) {
            if (file.exists())
                return true;
            else
                throw new UserException.CouldNotReadInputFile(file, "The interval file does not exist.");
        }

        if(file.exists())
            throw new UserException.CouldNotReadInputFile(file, String.format("The interval file %s does not have one of " +
                    "the supported extensions (.bed, .list, .picard, .interval_list, or .intervals). " +
                    "Please rename your file with the appropriate extension. If %s is NOT supposed to be a file, " +
                    "please move or rename the file at location %s", str, str, file.getAbsolutePath()));

        else return false;
    }

}


