package org.broadinstitute.sting.utils;

import edu.mit.broad.picard.directed.IntervalList;
import edu.mit.broad.picard.util.Interval;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jun 18, 2009
 * Time: 11:17:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenomeLocParser {
    private static Logger logger = Logger.getLogger(GenomeLocParser.class);

    private static final Pattern mPattern = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)?(\\+)?(-)?([\\d,]+)?$");  // matches case 3


    // --------------------------------------------------------------------------------------------------------------
    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    // --------------------------------------------------------------------------------------------------------------
    //public static Map<String, Integer> refContigOrdering = null;
    protected static SAMSequenceDictionary contigInfo = null;

    /**
     * do we have a contig ordering setup?
     *
     * @return true if the contig order is setup
     */
    public static boolean hasKnownContigOrdering() {
        return contigInfo != null;
    }

    /**
     * get the contig's SAMSequenceRecord
     *
     * @param contig the string name of the contig
     *
     * @return the sam sequence record
     */
    public static SAMSequenceRecord getContigInfo(final String contig) {
        return contigInfo.getSequence(contig);
    }

    /**
     * Returns the contig index of a specified string version of the contig
     *
     * @param contig the contig string
     *
     * @return the contig index, -1 if not found
     */
    public static int getContigIndex(final String contig) {
        if (contigInfo.getSequenceIndex(contig) == -1)
            Utils.scareUser(String.format("Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig));

        return contigInfo.getSequenceIndex(contig);
    }

    /**
     * set our internal reference contig order
     *
     * @param refFile the reference file
     *
     * @return true if we were successful
     */
    public static boolean setupRefContigOrdering(final ReferenceSequenceFile refFile) {
        return setupRefContigOrdering(refFile.getSequenceDictionary());
    }

    /**
     * setup our internal reference contig order
     *
     * @param seqDict the sequence dictionary
     *
     * @return true if we were successful
     */
    public static boolean setupRefContigOrdering(final SAMSequenceDictionary seqDict) {
        if (seqDict == null) { // we couldn't load the reference dictionary
            logger.info("Failed to load reference dictionary, falling back to lexicographic order for contigs");
            Utils.scareUser("Failed to load reference dictionary");
            return false;
        } else if (contigInfo == null) {
            contigInfo = seqDict;
            logger.debug(String.format("Prepared reference sequence contig dictionary"));
            for (SAMSequenceRecord contig : seqDict.getSequences()) {
                logger.debug(String.format(" %s (%d bp)", contig.getSequenceName(), contig.getSequenceLength()));
            }
        }
        return true;
    }

    /**
     * parse a genome location, from a location string
     *
     * @param str the string to parse
     *
     * @return a GenomeLoc representing the String
     */
    public static GenomeLoc parseGenomeLoc(final String str) {
        // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
        //System.out.printf("Parsing location '%s'%n", str);
        
        String contig = null;
        long start = 1;
        long stop = Integer.MAX_VALUE;
        boolean bad = false;

        Matcher match = mPattern.matcher(str);

        try {
            if (match.matches()) {
                contig = match.group(1);
                if (match.groupCount() == 5) {
                    start = parsePosition(match.group(2));
                    if (match.group(3) != null && match.group(3).equals("+") && match.group(5) == null) {
                        // do nothing
                    } else if (match.group(5) != null)
                        stop = parsePosition(match.group(5));
                    else if (match.group(5) == null && match.group(4) == null && match.group(3) == null)
                        stop = start;
                    else
                        bad = true;
                } else {
                    bad = true;
                }
            }
        }

        catch (Exception e) {
            bad = true;
        }

        if (bad)
            throw new StingException("Invalid Genome Location string: "); // + str);


        if (stop == Integer.MAX_VALUE && hasKnownContigOrdering())
            // lookup the actually stop position!
            stop = getContigInfo(contig).getSequenceLength();

        if (!isContigValid(contig))
            throw new MalformedGenomeLocException("Contig " + contig + " does not match any contig in the GATK sequence dictionary derived from the reference.");

        GenomeLoc loc = parseGenomeLoc(contig, start, stop);
        return loc;
    }

    /**
     * Useful utility function that parses a location string into a coordinate-order sorted
     * array of GenomeLoc objects
     *
     * @param str String representation of genome locs.  Null string corresponds to no filter.
     *
     * @return Array of GenomeLoc objects corresponding to the locations in the string, sorted by coordinate order
     */
    public static List<GenomeLoc> parseGenomeLocs(final String str) {
        // Null string means no filter.
        if (str == null) return null;

        // Of the form: loc1;loc2;...
        // Where each locN can be:
        // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
        try {
            List<GenomeLoc> locs = new ArrayList<GenomeLoc>();
            for (String loc : str.split(";"))
                locs.add(parseGenomeLoc(loc.trim()));
            Collections.sort(locs);
            //logger.info(String.format("Going to process %d locations", locs.length));
            locs = mergeOverlappingLocations(locs);
            logger.debug("Locations are:" + Utils.join(", ", locs));
            return locs;
        } catch (Exception e) { // TODO: fix this so that it passes the message from the exception, and doesn't print it out
            throw new StingException(String.format("Invalid locations string: %s, format is loc1;loc2; where each locN can be 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'", str), e);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Parsing string representations
    //
    // --------------------------------------------------------------------------------------------------------------
    private static long parsePosition(final String pos) {
        String x = pos.replaceAll(",", "");
        return Long.parseLong(x);
    }


    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     *
     * @return the list of merged locations
     */
    public static List<GenomeLoc> mergeOverlappingLocations(final List<GenomeLoc> raw) {
        logger.debug("  Raw locations are: " + Utils.join(", ", raw));
        if (raw.size() <= 1)
            return raw;
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<GenomeLoc>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while (it.hasNext()) {
                GenomeLoc curr = it.next();
                if (prev.contiguousP(curr)) {
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

    /**
     * Determines whether the given contig is valid with respect to the sequence dictionary
     * already installed in the GenomeLoc.
     *
     * @return True if the contig is valid.  False otherwise.
     */
    private static boolean isContigValid(String contig) {
        int contigIndex = contigInfo.getSequenceIndex(contig);
        return isSequenceIndexValid(contigIndex);
    }

    /**
     * Determines whether the given sequence index is valid with respect to the sequence dictionary.
     *
     * @param sequenceIndex sequence index
     *
     * @return True if the sequence index is valid, false otherwise.
     */
    private static boolean isSequenceIndexValid(int sequenceIndex) {
        return sequenceIndex >= 0 && sequenceIndex < contigInfo.size();

    }

    /**
     * Use this static constructor when the input data is under limited control (i.e. parsing user data).
     *
     * @param contig Contig to parse.
     * @param start  Starting point.
     * @param stop   Stop point.
     *
     * @return The genome location, or a MalformedGenomeLocException if unparseable.
     */
    public static GenomeLoc parseGenomeLoc(final String contig, long start, long stop) {
        if (!isContigValid(contig))
            throw new MalformedGenomeLocException("Contig " + contig + " does not match any contig in the GATK sequence dictionary derived from the reference.");
        return new GenomeLoc(contig, getContigIndex(contig), start, stop);
    }


    /**
     * Read a file of genome locations to process.
     * regions specified by the location string.  The string is of the form:
     * Of the form: loc1;loc2;...
     * Where each locN can be:
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
     *
     * @param file_name
     */
    public static List<GenomeLoc> intervalFileToList(final String file_name) {
        /**
         * first try to read it as an interval file since that's well structured
         * we'll fail quickly if it's not a valid file.  Then try to parse it as
         * a location string file
         */
        List<GenomeLoc> ret = null;
        try {
            File inputFile = new File(file_name);

            // sometimes we see an empty file passed as a parameter, if so return an empty list
            if (inputFile.exists() && inputFile.length() < 1) {
                if (GenomeAnalysisEngine.instance.getArguments().unsafe)
                    return new ArrayList<GenomeLoc>();
                else {
                    Utils.warnUser("The interval file " + file_name + " is empty. The GATK will continue processing but you " +
                            "may want to fix (or exclude) this file.");
                    return new ArrayList<GenomeLoc>();
                }
            }
            IntervalList il = IntervalList.fromFile(inputFile);

            // iterate through the list of merged intervals and add then as GenomeLocs
            ret = new ArrayList<GenomeLoc>();
            for (Interval interval : il.getUniqueIntervals()) {
                ret.add(new GenomeLoc(interval.getSequence(), getContigIndex(interval.getSequence()), interval.getStart(), interval.getEnd()));
            }
            return ret;

        } catch (Exception e) {
            try {
                xReadLines reader = new xReadLines(new File(file_name));
                List<String> lines = reader.readLines();
                reader.close();
                String locStr = Utils.join(";", lines);
                logger.debug("locStr: " + locStr);
                ret = parseGenomeLocs(locStr);
                return ret;
            } catch (Exception e2) {
                logger.error("Attempt to parse interval file in GATK format failed: " + e2.getMessage());
                e2.printStackTrace();
                throw new StingException("Unable to parse out interval file in either format", e);
            }
        }
    }

    /**
     * get the sequence name from a sequence index
     *
     * @param contigIndex get the contig index
     *
     * @return the string that represents that contig name
     */
    private static String getSequenceNameFromIndex(int contigIndex) {
        return GenomeLocParser.contigInfo.getSequence(contigIndex).getSequenceName();
    }

    /**
     * create a genome loc, given the contig name, start, and stop
     *
     * @param contig the contig name
     * @param start  the starting position
     * @param stop   the stop position
     *
     * @return a new genome loc
     */
    public static GenomeLoc createGenomeLoc(String contig, final long start, final long stop) {
        checkSetup();
        return verifyGenomeLoc(new GenomeLoc(contig, GenomeLocParser.getContigIndex(contig), start, stop));
    }

    /**
     * create a genome loc, given the contig index, start, and stop
     *
     * @param contigIndex the contig index
     * @param start       the start position
     * @param stop        the stop position
     *
     * @return a new genome loc
     */
    public static GenomeLoc createGenomeLoc(int contigIndex, final long start, final long stop) {
        checkSetup();
        if (start < 0) {
            throw new StingException("Bad start position " + start);
        }
        if (stop < -1) {
            throw new StingException("Bad stop position " + stop);
        }    // a negative -1 indicates it's not a meaningful end position


        return new GenomeLoc(getSequenceNameFromIndex(contigIndex), contigIndex, start, stop);
    }

    /**
     * create a genome loc, given a read
     *
     * @param read
     *
     * @return
     */
    public static GenomeLoc createGenomeLoc(final SAMRecord read) {
        checkSetup();
        return verifyGenomeLoc(new GenomeLoc(read.getReferenceName(), read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd()));
    }


    /**
     * create a new genome loc, given the contig position, and a single position
     *
     * @param contig the contig name
     * @param pos    the postion
     *
     * @return a genome loc representing a single base at the specified postion on the contig
     */
    public static GenomeLoc createGenomeLoc(final int contig, final long pos) {
        checkSetup();
        return verifyGenomeLoc(new GenomeLoc(getSequenceNameFromIndex(contig), contig, pos, pos));
    }

    /**
     * create a new genome loc, given the contig name, and a single position
     *
     * @param contig the contig name
     * @param pos    the postion
     *
     * @return a genome loc representing a single base at the specified postion on the contig
     */
    public static GenomeLoc createGenomeLoc(final String contig, final long pos) {
        checkSetup();
        return verifyGenomeLoc(new GenomeLoc(contig, GenomeLocParser.getContigIndex(contig), pos, pos));
    }

    public static GenomeLoc createGenomeLoc(final GenomeLoc toCopy) {
        checkSetup();
        return verifyGenomeLoc(new GenomeLoc(toCopy.getContig(), toCopy.getContigIndex(), toCopy.getStart(), toCopy.getStop()));
    }

    /**
     * verify the specified genome loc is valid, if it's not, throw an exception
     *
     * @param toReturn teh genome loc we're about to return
     *
     * @return the genome loc if it's valid, otherwise we throw an exception
     */
    private static GenomeLoc verifyGenomeLoc(GenomeLoc toReturn) {
        if (toReturn.getStart() < 0) {
            throw new StingException("Parameters to GenomeLocParser are incorrect: the start position is less than 0");
        }
        if ((toReturn.getStop() != -1) && (toReturn.getStop() < 0)) {
            throw new StingException("Parameters to GenomeLocParser are incorrect: the stop position is less than 0");
        }
        if (toReturn.getContigIndex() < 0) {
            throw new StingException("Parameters to GenomeLocParser are incorrect: the contig index is less than 0");
        }
        if (toReturn.getContigIndex() >= contigInfo.getSequences().size()) {
            throw new StingException("Parameters to GenomeLocParser are incorrect: the contig index is greater then the stored sequence count");

        }
        return toReturn;

    }


    /**
     * Move this Genome loc to the next contig, with a start
     * and stop of 1.
     *
     * @return true if we are not out of contigs, otherwise false if we're
     *         at the end of the genome (no more contigs to jump to).
     */
    public static GenomeLoc toNextContig(GenomeLoc current) {
        if (current.getContigIndex() + 1 >= contigInfo.getSequences().size()) {
            return null;
        } else
            return verifyGenomeLoc(new GenomeLoc(getSequenceNameFromIndex(current.getContigIndex() + 1), current.getContigIndex() + 1, 1, 1));
    }

    /**
     * create a new genome loc, given an old location and a new contig
     *
     * @param loc    the old location
     * @param contig the new contig to set
     *
     * @return a new genome loc with an updated contig name and index
     */
    public static GenomeLoc setContig(GenomeLoc loc, String contig) {
        checkSetup();

        int index = -1;
        if ((index = contigInfo.getSequenceIndex(contig)) < 0) {
            throw new StingException("Contig name ( " + contig + " ) not in the set sequence dictionary.");
        }
        return verifyGenomeLoc(new GenomeLoc(contig, index, loc.start, loc.getStop()));
    }

    /**
     * Sets contig index. UNSAFE since it 1) does NOT update contig name; 2) does not validate the index
     *
     * @param contig
     */
    public static GenomeLoc setContigIndex(GenomeLoc loc, int contig) {
        checkSetup();
        if ((contig >= GenomeLocParser.contigInfo.getSequences().size()) || (contig < 0)) {
            throw new StingException("Contig index ( " + contig + " ) is not in the sequence dictionary set.");
        }
        return verifyGenomeLoc(new GenomeLoc(GenomeLocParser.contigInfo.getSequence(contig).getSequenceName(), contig, loc.start, loc.getStop()));
    }


    /**
     * create a new genome loc from an existing loc, with a new start position
     *
     * @param loc   the old location
     * @param start a new start position
     *
     * @return the newly created genome loc
     */
    public static GenomeLoc setStart(GenomeLoc loc, long start) {
        checkSetup();
        if (loc.getContigIndex() < 0 || loc.getContigIndex() >= contigInfo.getSequences().size()) {
            throw new StingException("Genome loc passed in to setStart has a contig index outside the range of our current sequence dictionary");
        }
        int length = GenomeLocParser.contigInfo.getSequences().get(loc.getContigIndex()).getSequenceLength();
        if (start > length) {
            throw new StingException("start value of " + start + " is greater than the contig length, and is not -1. (length = " + length + ").");
        }
        return verifyGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), start, loc.getStop()));
    }

    /**
     * create a new genome loc from an existing loc, with a new stop position
     *
     * @param loc  the old location
     * @param stop a new stop position
     *
     * @return
     */
    public static GenomeLoc setStop(GenomeLoc loc, long stop) {
        checkSetup();
        if (loc.getContigIndex() < 0 || loc.getContigIndex() >= contigInfo.getSequences().size()) {
            throw new StingException("Genome loc passed in to setStart has a contig index outside the range of our current sequence dictionary");
        }
        if ((stop != -1) && (stop > GenomeLocParser.contigInfo.getSequences().get(loc.getContigIndex()).getSequenceLength())) {
            throw new StingException("stop value of " + stop + " is greater than the contig length, and is not -1.");
        }
        return verifyGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start, stop));
    }

    /**
     * return a new genome loc, with an incremented position
     *
     * @param loc the old location
     *
     * @return a new genome loc
     */
    public static GenomeLoc incPos(GenomeLoc loc) {
        return incPos(loc, 1);
    }

    /**
     * return a new genome loc, with an incremented position
     *
     * @param loc the old location
     * @param by  how much to move the start and stop by
     *
     * @return a new genome loc
     */
    public static GenomeLoc incPos(GenomeLoc loc, long by) {
        return verifyGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start + by, loc.stop + by));
    }

    /**
     * create a new genome loc with an incremented position
     *
     * @param loc the location
     *
     * @return a new genome loc
     */
    public static GenomeLoc nextLoc(GenomeLoc loc) {
        return incPos(loc);
    }

    /** check to make sure that we've setup the contig information */
    private static void checkSetup() {
        if (contigInfo == null) {
            throw new StingException("The GenomeLocParser hasn't been setup with a contig sequence yet");
        }
    }

    /**
     * compare two contig names, in the current context
     *
     * @param firstContig
     * @param secondContig
     *
     * @return
     */
    public static int compareContigs(String firstContig, String secondContig) {
        checkSetup();
        Integer ref1 = GenomeLocParser.getContigIndex(firstContig);
        Integer ref2 = GenomeLocParser.getContigIndex(secondContig);
        return ref1.compareTo(ref2);

    }
}
