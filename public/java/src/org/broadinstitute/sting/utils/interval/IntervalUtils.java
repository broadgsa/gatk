package org.broadinstitute.sting.utils.interval;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
     * @return an unsorted, unmerged representation of the given intervals.  Null is used to indicate that all intervals should be used.
     */
    public static List<GenomeLoc> parseIntervalArguments(GenomeLocParser parser, List<String> argList) {
        List<GenomeLoc> rawIntervals = new ArrayList<GenomeLoc>();    // running list of raw GenomeLocs

        if (argList != null) { // now that we can be in this function if only the ROD-to-Intervals was provided, we need to
                               // ensure that the arg list isn't null before looping.
            for (String argument : argList) {
                rawIntervals.addAll(parseIntervalArguments(parser, argument));
            }
        }

        return rawIntervals;
    }

    public static List<GenomeLoc> parseIntervalArguments(GenomeLocParser parser, String arg) {
        List<GenomeLoc> rawIntervals = new ArrayList<GenomeLoc>();    // running list of raw GenomeLocs

        // separate argument on semicolon first
        for (String fileOrInterval : arg.split(";")) {
            // if any argument is 'unmapped', "parse" it to a null entry.  A null in this case means 'all the intervals with no alignment data'.
            if (isUnmapped(fileOrInterval))
                rawIntervals.add(GenomeLoc.UNMAPPED);
            // if it's a file, add items to raw interval list
            else if (isIntervalFile(fileOrInterval)) {
                try {
                    rawIntervals.addAll(intervalFileToList(parser, fileOrInterval));
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

        return rawIntervals;
    }

    /**
     * Read a file of genome locations to process. The file may be in BED, Picard,
     * or GATK interval format.
     *
     * @param glParser   GenomeLocParser
     * @param file_name  interval file
     * @return List<GenomeLoc> List of Genome Locs that have been parsed from file
     */
    public static List<GenomeLoc> intervalFileToList(final GenomeLocParser glParser, final String file_name) {
        // try to open file
        File inputFile = new File(file_name);
        List<GenomeLoc> ret = new ArrayList<GenomeLoc>();

        // case: BED file
        if ( file_name.toUpperCase().endsWith(".BED") ) {
            // this is now supported in Tribble
            throw new ReviewedStingException("BED files must be parsed through Tribble; parsing them as intervals through the GATK engine is no longer supported");
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
                throw new UserException.BadInput("The INTERSECTION of your -L options produced no intervals.");

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
     * computes whether the test interval list is equivalent to master.  To be equivalent, test must
     * contain GenomeLocs covering every base in master, exactly once.  Note that this algorithm
     * assumes that master genomelocs are all discontiguous (i.e., we don't have locs like 1-3 and 4-6 but
     * rather just 1-6).  In order to use this algorithm with contiguous genomelocs first merge them.  The algorithm
     * doesn't assume that test has discontinuous genomelocs.
     *
     * Returns a null string if there are no differences, otherwise returns a string describing the difference
     * (useful for UnitTests).  Assumes both lists are sorted
     */
    public static final String equateIntervals(List<GenomeLoc> masterArg, List<GenomeLoc> testArg) {
        LinkedList<GenomeLoc> master = new LinkedList<GenomeLoc>(masterArg);
        LinkedList<GenomeLoc> test = new LinkedList<GenomeLoc>(testArg);

        while ( ! master.isEmpty() ) { // there's still unchecked bases in master
            final GenomeLoc masterHead = master.pop();
            final GenomeLoc testHead = test.pop();

            if ( testHead.overlapsP(masterHead) ) {
                // remove the parts of test that overlap master, and push the remaining
                // parts onto master for further comparison.
                for ( final GenomeLoc masterPart : Utils.reverse(masterHead.subtract(testHead)) ) {
                    master.push(masterPart);
                }
            } else {
                // testHead is incompatible with masterHead, so we must have extra bases in testHead
                // that aren't in master
                return "Incompatible locs detected masterHead=" + masterHead + ", testHead=" + testHead;
            }
        }

        if ( test.isEmpty() ) // everything is equal
            return null; // no differences
        else
            return "Remaining elements found in test: first=" + test.peek();
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
     * Splits an interval list into multiple sublists.
     * @param locs The genome locs to split.
     * @param splits The stop points for the genome locs returned by splitFixedIntervals.
     * @return A list of lists of genome locs, split according to splits
     */
    public static List<List<GenomeLoc>> splitIntervalsToSubLists(List<GenomeLoc> locs, List<Integer> splits) {
        int locIndex = 1;
        int start = 0;
        List<List<GenomeLoc>> sublists = new ArrayList<List<GenomeLoc>>(splits.size());
        for (Integer stop: splits) {
            List<GenomeLoc> curList = new ArrayList<GenomeLoc>();
            for (int i = start; i < stop; i++)
                curList.add(locs.get(i));
            start = stop;
            sublists.add(curList);
        }

        return sublists;
    }


    /**
     * Splits an interval list into multiple files.
     * @param fileHeader The sam file header.
     * @param splits Pre-divided genome locs returned by splitFixedIntervals.
     * @param scatterParts The output interval lists to write to.
     */
    public static void scatterFixedIntervals(SAMFileHeader fileHeader, List<List<GenomeLoc>> splits, List<File> scatterParts) {
        if (splits.size() != scatterParts.size())
            throw new UserException.BadArgumentValue("splits", String.format("Split points %d does not equal the number of scatter parts %d.", splits.size(), scatterParts.size()));

        int fileIndex = 0;
        int locIndex = 1;
        for (final List<GenomeLoc> split : splits) {
            IntervalList intervalList = new IntervalList(fileHeader);
            for (final GenomeLoc loc : split)
                intervalList.add(toInterval(loc, locIndex++));
            intervalList.write(scatterParts.get(fileIndex++));
        }
    }

    /**
     * Splits the genome locs up by size.
     * @param locs Genome locs to split.
     * @param numParts Number of parts to split the locs into.
     * @return The stop points to split the genome locs.
     */
    public static List<List<GenomeLoc>> splitFixedIntervals(List<GenomeLoc> locs, int numParts) {
        if (locs.size() < numParts)
            throw new UserException.BadArgumentValue("scatterParts", String.format("Cannot scatter %d locs into %d parts.", locs.size(), numParts));
        final long locsSize = intervalSize(locs);
        final List<Integer> splitPoints = new ArrayList<Integer>();
        addFixedSplit(splitPoints, locs, locsSize, 0, locs.size(), numParts);
        Collections.sort(splitPoints);
        splitPoints.add(locs.size());
        return splitIntervalsToSubLists(locs, splitPoints);
    }

    @Requires({"locs != null", "numParts > 0"})
    @Ensures("result != null")
    public static List<List<GenomeLoc>> splitLocusIntervals(List<GenomeLoc> locs, int numParts) {
        // the ideal size of each split
        final long bp = IntervalUtils.intervalSize(locs);
        final long idealSplitSize = Math.max((long)Math.floor(bp / (1.0*numParts)), 1);

        // algorithm:
        // split = ()
        // set size = 0
        // pop the head H off locs.
        // If size + size(H) < splitSize:
        //      add H to split, continue
        // If size + size(H) == splitSize:
        //      done with split, put in splits, restart
        // if size + size(H) > splitSize:
        //      cut H into two pieces, first of which has splitSize - size bp
        //      push both pieces onto locs, continue
        // The last split is special -- when you have only one split left, it gets all of the remaining locs
        // to deal with rounding issues
        final List<List<GenomeLoc>> splits = new ArrayList<List<GenomeLoc>>(numParts);

        LinkedList<GenomeLoc> locsLinkedList = new LinkedList<GenomeLoc>(locs);
        while ( ! locsLinkedList.isEmpty() ) {
            if ( splits.size() + 1 == numParts ) {
                // the last one gets all of the remaining parts
                splits.add(new ArrayList<GenomeLoc>(locsLinkedList));
                locsLinkedList.clear();
            } else {
                final SplitLocusRecursive one = splitLocusIntervals1(locsLinkedList, idealSplitSize);
                splits.add(one.split);
                locsLinkedList = one.remaining;
            }
        }

        return splits;
    }

    @Requires({"remaining != null", "!remaining.isEmpty()", "idealSplitSize > 0"})
    @Ensures({"result != null"})
    final static SplitLocusRecursive splitLocusIntervals1(LinkedList<GenomeLoc> remaining, long idealSplitSize) {
        final List<GenomeLoc> split = new ArrayList<GenomeLoc>();
        long size = 0;

        while ( ! remaining.isEmpty() ) {
            GenomeLoc head = remaining.pop();
            final long newSize = size + head.size();

            if ( newSize == idealSplitSize ) {
                split.add(head);
                break; // we are done
            } else if ( newSize > idealSplitSize ) {
                final long remainingBp = idealSplitSize - size;
                final long cutPoint = head.getStart() + remainingBp;
                GenomeLoc[] parts = head.split((int)cutPoint);
                remaining.push(parts[1]);
                remaining.push(parts[0]);
                // when we go around, head.size' = idealSplitSize - size
                // so newSize' = splitSize + head.size' = size + (idealSplitSize - size) = idealSplitSize
            } else {
                split.add(head);
                size = newSize;
            }
        }

        return new SplitLocusRecursive(split, remaining);
    }

    private final static class SplitLocusRecursive {
        final List<GenomeLoc> split;
        final LinkedList<GenomeLoc> remaining;

        @Requires({"split != null", "remaining != null"})
        private SplitLocusRecursive(final List<GenomeLoc> split, final LinkedList<GenomeLoc> remaining) {
            this.split = split;
            this.remaining = remaining;
        }
    }

    public static List<GenomeLoc> flattenSplitIntervals(List<List<GenomeLoc>> splits) {
        final List<GenomeLoc> locs = new ArrayList<GenomeLoc>();
        for ( final List<GenomeLoc> split : splits )
            locs.addAll(split);
        return locs;
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

    public static final long intervalSize(final List<GenomeLoc> locs) {
        long size = 0;
        for ( final GenomeLoc loc : locs )
            size += loc.size();
        return size;
    }
}
