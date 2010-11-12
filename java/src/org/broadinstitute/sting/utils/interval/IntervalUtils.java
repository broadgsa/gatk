package org.broadinstitute.sting.utils.interval;

import net.sf.picard.util.IntervalList;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceDataSource;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
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
     * @param parser Genome loc parser.
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
     * @param parser Genome loc parser for the intervals.
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

    /**
     * Returns the list of GenomeLocs from the list of intervals.
     * @param referenceSource The reference for the intervals.
     * @param intervals The interval as strings or file paths.
     * @return The list of GenomeLocs.
     */
    private static List<GenomeLoc> parseIntervalArguments(ReferenceDataSource referenceSource, List<String> intervals) {
        GenomeLocParser parser = new GenomeLocParser(referenceSource.getReference());
        GenomeLocSortedSet locs;
        // TODO: Abstract genome analysis engine has richer logic for parsing.  We need to use it!
        if (intervals.size() == 0) {
            locs = GenomeLocSortedSet.createSetFromSequenceDictionary(referenceSource.getReference().getSequenceDictionary());
        } else {
            locs = new GenomeLocSortedSet(parser, IntervalUtils.parseIntervalArguments(parser, intervals, false));
        }
        if (locs == null || locs.size() == 0)
            throw new UserException.MalformedFile("Intervals are empty: " + Utils.join(", ", intervals));
        return locs.toList();
    }

    /**
     * Returns the list of contigs from the list of intervals.
     * @param reference The reference for the intervals.
     * @return The list of contig names.
     */
    public static List<String> distinctContigs(File reference) {
        return distinctContigs(reference, Collections.<String>emptyList());
    }

    /**
     * Returns the list of contigs from the list of intervals.
     * @param reference The reference for the intervals.
     * @param intervals The interval as strings or file paths.
     * @return The list of contig names.
     */
    public static List<String> distinctContigs(File reference, List<String> intervals) {
        ReferenceDataSource referenceSource = new ReferenceDataSource(reference);
        List<GenomeLoc> locs = parseIntervalArguments(referenceSource, intervals);
        String contig = null;
        List<String> contigs = new ArrayList<String>();
        for (GenomeLoc loc: locs) {
            if (contig == null || !contig.equals(loc.getContig())) {
                contig = loc.getContig();
                contigs.add(contig);
            }
        }
        return contigs;
    }

    /**
     * Splits an interval list into multiple files.
     * @param reference The reference for the intervals.
     * @param intervals The interval as strings or file paths.
     * @param scatterParts The output interval lists to write to.
     * @param splitByContig If true then one contig will not be written to multiple files.
     */
    public static void scatterIntervalArguments(File reference, List<String> intervals, List<File> scatterParts, boolean splitByContig) {
        ReferenceDataSource referenceSource = new ReferenceDataSource(reference);
        List<GenomeLoc> locs = parseIntervalArguments(referenceSource, intervals);
        SAMFileHeader fileHeader = new SAMFileHeader();
        fileHeader.setSequenceDictionary(referenceSource.getReference().getSequenceDictionary());

        IntervalList intervalList = null;
        int fileIndex = -1;
        int locIndex = 0;

        if (splitByContig) {
            String contig = null;
            for (GenomeLoc loc: locs) {
                // If there are still more files to write and the contig doesn't match...
                if ((fileIndex+1 < scatterParts.size()) && (contig == null || !contig.equals(loc.getContig()))) {
                    // Then close the current file and start a new one.
                    if (intervalList != null) {
                        intervalList.write(scatterParts.get(fileIndex));
                        intervalList = null;
                    }
                    fileIndex++;
                    contig = loc.getContig();
                }
                if (intervalList == null)
                    intervalList = new IntervalList(fileHeader);
                intervalList.add(toInterval(loc, ++locIndex));
            }
            if (intervalList != null)
                intervalList.write(scatterParts.get(fileIndex));
        } else {
            int locsPerFile = locs.size() / scatterParts.size();
            int locRemainder = locs.size() % scatterParts.size();

            // At the start, put an extra loc per file
            locsPerFile++;
            int locsLeftFile = 0;

            for (GenomeLoc loc: locs) {
                if (locsLeftFile == 0) {
                    if (intervalList != null)
                        intervalList.write(scatterParts.get(fileIndex));

                    fileIndex++;
                    intervalList = new IntervalList(fileHeader);

                    // When we have put enough locs into each file,
                    // reduce the number of locs per file back
                    // to the original calculated value.
                    if (fileIndex == locRemainder)
                        locsPerFile -= 1;
                    locsLeftFile = locsPerFile;
                }
                locsLeftFile -= 1;
                intervalList.add(toInterval(loc, ++locIndex));
            }
            if (intervalList != null)
                intervalList.write(scatterParts.get(fileIndex));
        }
        if ((fileIndex + 1) != scatterParts.size())
            throw new UserException.BadArgumentValue("scatterParts", String.format("Only able to write contigs into %d of %d files.", fileIndex + 1, scatterParts.size()));
    }

    /**
     * Converts a GenomeLoc to a picard interval.
     * @param loc The GenomeLoc.
     * @param locIndex The loc index for use in the file.
     * @return The picard interval.
     */
    private static net.sf.picard.util.Interval toInterval(GenomeLoc loc, int locIndex) {
        return new net.sf.picard.util.Interval(loc.getContig(), loc.getStart(), loc.getStop(), false, "interval_" + locIndex);
    }
}
