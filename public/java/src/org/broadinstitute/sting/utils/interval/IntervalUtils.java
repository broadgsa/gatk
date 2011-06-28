package org.broadinstitute.sting.utils.interval;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.bed.BedParser;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.IOException;
import java.util.*;
import java.io.File;

/**
 * Parse text representations of interval strings that
 * can appear in Sting-based applications.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalUtils {
    private static Logger logger = Logger.getLogger(IntervalUtils.class);

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

                // separate argument on semicolon first
                for (String fileOrInterval : argument.split(";")) {
                    // if any interval argument is '-L all', consider all loci by returning no intervals
                    if (fileOrInterval.trim().toLowerCase().equals("all")) {
                        if (argList.size() != 1) {
                            // throw error if '-L all' is not only interval - potentially conflicting commands
                            throw new UserException.CommandLineException(String.format("Conflicting arguments: Intervals given along with \"-L all\""));
                        }
                        return null;
                    }
                    // if any argument is 'unmapped', "parse" it to a null entry.  A null in this case means 'all the intervals with no alignment data'.
                    else if (isUnmapped(fileOrInterval))
                        rawIntervals.add(GenomeLoc.UNMAPPED);
                    // if it's a file, add items to raw interval list
                    else if (isIntervalFile(fileOrInterval)) {
                        try {
                            rawIntervals.addAll(intervalFileToList(parser, fileOrInterval, allowEmptyIntervalList));
                        }
                        catch ( UserException.MalformedGenomeLoc e ) {
                            throw e;
                        }
                        catch ( Exception e ) {
                            throw new UserException.MalformedFile(fileOrInterval, "Interval file could not be parsed in any supported format.", e);
                        }
                    }

                        // otherwise treat as an interval -> parse and add to raw interval list
                    else {
                        rawIntervals.add(parser.parseGenomeLoc(fileOrInterval));
                    }
                }
            }
        }

        return rawIntervals;
    }

        /**
     * Read a file of genome locations to process. The file may be in BED, Picard,
     * or GATK interval format.
     *
     * @param file_name interval file
     * @param allowEmptyIntervalList if false an exception will be thrown for files that contain no intervals
     * @return List<GenomeLoc> List of Genome Locs that have been parsed from file
     */
    public static List<GenomeLoc> intervalFileToList(final GenomeLocParser glParser, final String file_name, boolean allowEmptyIntervalList) {
        // try to open file
        File inputFile = new File(file_name);
        List<GenomeLoc> ret = new ArrayList<GenomeLoc>();

        // case: BED file
        if (file_name.toUpperCase().endsWith(".BED")) {
            BedParser parser = new BedParser(glParser,inputFile);
            ret.addAll(parser.getLocations());
        }
        else {
            /**
             * IF not a BED file:
             * first try to read it as a Picard interval file since that's well structured
             * we'll fail quickly if it's not a valid file.
             */
            boolean isPicardInterval = false;
            try {
                // Note: Picard will skip over intervals with contigs not in the sequence dictionary
                IntervalList il = IntervalList.fromFile(inputFile);
                isPicardInterval = true;

                int nInvalidIntervals = 0;
                for (Interval interval : il.getIntervals()) {
                    if ( glParser.isValidGenomeLoc(interval.getSequence(), interval.getStart(), interval.getEnd(), true))
                        ret.add(glParser.createGenomeLoc(interval.getSequence(), interval.getStart(), interval.getEnd(), true));
                    else {
                        nInvalidIntervals++;
                    }
                }
                if ( nInvalidIntervals > 0 )
                    logger.warn("Ignoring " + nInvalidIntervals + " invalid intervals from " + inputFile);
            }

            // if that didn't work, try parsing file as a GATK interval file
            catch (Exception e) {
                if ( isPicardInterval ) // definitely a picard file, but we failed to parse
                    throw new UserException.CouldNotReadInputFile(inputFile, e);
                else {
                    try {
                        XReadLines reader = new XReadLines(new File(file_name));
                        for(String line: reader) {
                            if ( line.trim().length() > 0 ) {
                                ret.add(glParser.parseGenomeLoc(line));
                            }
                        }
                        reader.close();
                    }
                    catch (IOException e2) {
                        throw new UserException.CouldNotReadInputFile(inputFile, e2);
                    }
                }
            }
        }

        if ( ret.isEmpty() && ! allowEmptyIntervalList ) {
            throw new UserException("The interval file " + inputFile.getAbsolutePath() + " contains no intervals " +
                                    "that could be parsed, and the unsafe operation ALLOW_EMPTY_INTERVAL_LIST has " +
                                    "not been enabled");
        }

        return ret;
    }

    /**
     * Returns true if the interval string is the "unmapped" interval
     * @param interval Interval to check
     * @return true if the interval string is the "unmapped" interval
     */
    public static boolean isUnmapped(String interval) {
        return (interval != null && interval.trim().toLowerCase().equals("unmapped"));
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

        //if we have an empty list, throw an exception.  If they specified intersection and there are no items, this is bad.
        if (retList.size() == 0)
                throw new UserException.BadInput("The INTERSECTION of your -BTI and -L options produced no intervals.");

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
        intervals = mergeIntervalLocations(intervals, mergingRule);

        return GenomeLocSortedSet.createSetFromList(parser,intervals);
    }

    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(String str) {
        return isIntervalFile(str, true);
    }

    /**
     * Check if string argument was intented as a file
     * Accepted file extensions: .bed .list, .picard, .interval_list, .intervals.
     * @param str token to identify as a filename.
     * @param checkExists if true throws an exception if the file doesn't exist.
     * @return true if the token looks like a filename, or false otherwise.
     */
    public static boolean isIntervalFile(String str, boolean checkExists) {
        // should we define list of file extensions as a public array somewhere?
        // is regex or endsiwth better?
        File file = new File(str);
        if (str.toUpperCase().endsWith(".BED") || str.toUpperCase().endsWith(".LIST") ||
                str.toUpperCase().endsWith(".PICARD") || str.toUpperCase().endsWith(".INTERVAL_LIST")
                || str.toUpperCase().endsWith(".INTERVALS")) {
            if (!checkExists)
                return true;
            else if (file.exists())
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
     * Returns a map of contig names with their sizes.
     * @param reference The reference for the intervals.
     * @return A map of contig names with their sizes.
     */
    public static Map<String, Long> getContigSizes(File reference) {
        ReferenceDataSource referenceSource = new ReferenceDataSource(reference);
        List<GenomeLoc> locs = GenomeLocSortedSet.createSetFromSequenceDictionary(referenceSource.getReference().getSequenceDictionary()).toList();
        Map<String, Long> lengths = new LinkedHashMap<String, Long>();
        for (GenomeLoc loc: locs)
            lengths.put(loc.getContig(), loc.size());
        return lengths;
    }

    /**
     * Counts the number of interval files an interval list can be split into using scatterIntervalArguments.
     * @param locs The genome locs.
     * @return The maximum number of parts the intervals can be split into.
     */
    public static int countContigIntervals(List<GenomeLoc> locs) {
        int maxFiles = 0;
        String contig = null;
        for (GenomeLoc loc: locs) {
            if (contig == null || !contig.equals(loc.getContig())) {
                maxFiles++;
                contig = loc.getContig();
            }
        }
        return maxFiles;
    }

    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param locs The genome locs to split.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterContigIntervals(SAMFileHeader fileHeader, List<GenomeLoc> locs, List<File> scatterParts) {
        IntervalList intervalList = null;
        int fileIndex = -1;
        int locIndex = 0;
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
        if ((fileIndex + 1) != scatterParts.size())
            throw new UserException.BadArgumentValue("scatterParts", String.format("Only able to write contigs into %d of %d files.", fileIndex + 1, scatterParts.size()));
    }

    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param locs The genome locs to split.
     * @param splits The stop points for the genome locs returned by splitFixedIntervals.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterFixedIntervals(SAMFileHeader fileHeader, List<GenomeLoc> locs, List<Integer> splits, List<File> scatterParts) {
        if (splits.size() != scatterParts.size())
            throw new UserException.BadArgumentValue("splits", String.format("Split points %d does not equal the number of scatter parts %d.", splits.size(), scatterParts.size()));
        int fileIndex = 0;
        int locIndex = 1;
        int start = 0;
        for (Integer stop: splits) {
            IntervalList intervalList = new IntervalList(fileHeader);
            for (int i = start; i < stop; i++)
                intervalList.add(toInterval(locs.get(i), locIndex++));
            intervalList.write(scatterParts.get(fileIndex++));
            start = stop;
        }
    }

    /**
     * Splits the genome locs up by size.
     * @param locs Genome locs to split.
     * @param numParts Number of parts to split the locs into.
     * @return The stop points to split the genome locs.
     */
    public static List<Integer> splitFixedIntervals(List<GenomeLoc> locs, int numParts) {
        if (locs.size() < numParts)
            throw new UserException.BadArgumentValue("scatterParts", String.format("Cannot scatter %d locs into %d parts.", locs.size(), numParts));
        long locsSize = 0;
        for (GenomeLoc loc: locs)
            locsSize += loc.size();
        List<Integer> splitPoints = new ArrayList<Integer>();
        addFixedSplit(splitPoints, locs, locsSize, 0, locs.size(), numParts);
        Collections.sort(splitPoints);
        splitPoints.add(locs.size());
        return splitPoints;
    }

    private static void addFixedSplit(List<Integer> splitPoints, List<GenomeLoc> locs, long locsSize, int startIndex, int stopIndex, int numParts) {
        if (numParts < 2)
            return;
        int halfParts = (numParts + 1) / 2;
        Pair<Integer, Long> splitPoint = getFixedSplit(locs, locsSize, startIndex, stopIndex, halfParts, numParts - halfParts);
        int splitIndex = splitPoint.first;
        long splitSize = splitPoint.second;
        splitPoints.add(splitIndex);
        addFixedSplit(splitPoints, locs, splitSize, startIndex, splitIndex, halfParts);
        addFixedSplit(splitPoints, locs, locsSize - splitSize, splitIndex, stopIndex, numParts - halfParts);
    }

    private static Pair<Integer, Long> getFixedSplit(List<GenomeLoc> locs, long locsSize, int startIndex, int stopIndex, int minLocs, int maxLocs) {
        int splitIndex = startIndex;
        long splitSize = 0;
        for (int i = 0; i < minLocs; i++) {
            splitSize += locs.get(splitIndex).size();
            splitIndex++;
        }
        long halfSize = locsSize / 2;
        while (splitIndex < (stopIndex - maxLocs) && splitSize < halfSize) {
            splitSize += locs.get(splitIndex).size();
            splitIndex++;
        }
        return new Pair<Integer, Long>(splitIndex, splitSize);
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

    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     * @param rule the merging rule we're using
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeIntervalLocations(final List<GenomeLoc> raw, IntervalMergingRule rule) {
        if (raw.size() <= 1)
            return raw;
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<GenomeLoc>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                GenomeLoc curr = it.next();
                if (prev.overlapsP(curr)) {
                    prev = prev.merge(curr);
                } else if (prev.contiguousP(curr) && rule == IntervalMergingRule.ALL) {
                    prev = prev.merge(curr);
                } else {
                    merged.add(prev);
                    prev = curr;
                }
            }
            merged.add(prev);
            return merged;
        }
    }
}
