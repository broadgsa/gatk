/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.utils.bed.BedParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.text.XReadLines;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jun 18, 2009
 * Time: 11:17:01 PM
 * To change this template use File | Settings | File Templates.
 */
public class GenomeLocParser {
    private static Logger logger = Logger.getLogger(GenomeLocParser.class);

    // --------------------------------------------------------------------------------------------------------------
    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    // --------------------------------------------------------------------------------------------------------------
    //public static Map<String, Integer> refContigOrdering = null;
    protected SAMSequenceDictionary contigInfo = null;

    /**
     * set our internal reference contig order
     * @param refFile the reference file
     */
    public GenomeLocParser(final ReferenceSequenceFile refFile) {
        this(refFile.getSequenceDictionary());
    }

    public GenomeLocParser(SAMSequenceDictionary seqDict) {
        if (seqDict == null) { // we couldn't load the reference dictionary
            //logger.info("Failed to load reference dictionary, falling back to lexicographic order for contigs");
            throw new UserException.CommandLineException("Failed to load reference dictionary");
        } else if (contigInfo == null) {
            contigInfo = seqDict;
            logger.debug(String.format("Prepared reference sequence contig dictionary"));
            for (SAMSequenceRecord contig : seqDict.getSequences()) {
                logger.debug(String.format(" %s (%d bp)", contig.getSequenceName(), contig.getSequenceLength()));
            }
        }
    }

    /**
     * get the contig's SAMSequenceRecord
     *
     * @param contig the string name of the contig
     *
     * @return the sam sequence record
     */
    public SAMSequenceRecord getContigInfo(final String contig) {
        return contigInfo.getSequence(contig);
    }

    /**
     * Returns the contig index of a specified string version of the contig
     *
     * @param contig the contig string
     * @param exceptionOut in some cases we don't want to exception out if the contig isn't valid
     *
     * @return the contig index, -1 if not found
     */
    public int getContigIndex(final String contig, boolean exceptionOut) {
        if (contigInfo.getSequenceIndex(contig) == -1 && exceptionOut)
            throw new UserException.CommandLineException(String.format("Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig));

        return contigInfo.getSequenceIndex(contig);
    }

    /**
     * parse a genome interval, from a location string
     *
     * Performs interval-style validation:
     *
     * contig is valid; start and stop less than the end; start <= sto
     * @param str the string to parse
     *
     * @return a GenomeLoc representing the String
     *
     */

    public GenomeLoc parseGenomeInterval(final String str) {
         GenomeLoc ret = parseGenomeLoc(str);
         exceptionOnInvalidGenomeLocBounds(ret);
         return ret;
     }

    /**
     * parse a genome location, from a location string
     *
     * Performs read-style validation:
     * checks that start and stop are positive, start < stop, and the contig is valid
     * does not check that genomeLoc is actually on the contig
     *
     * @param str the string to parse
     *
     * @return a GenomeLoc representing the String
     *
     */
    public GenomeLoc parseGenomeLoc(final String str) {
        // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
        //System.out.printf("Parsing location '%s'%n", str);

        String contig = null;
        int start = 1;
        int stop = -1;

        final int colonIndex = str.indexOf(":");
        if(colonIndex == -1) {
            contig = str.substring(0, str.length());  // chr1
            stop = Integer.MAX_VALUE;
        } else {
            contig = str.substring(0, colonIndex);
            final int dashIndex = str.indexOf('-', colonIndex);
            try {
                if(dashIndex == -1) {
                    if(str.charAt(str.length() - 1) == '+') {
                        start = parsePosition(str.substring(colonIndex + 1, str.length() - 1));  // chr:1+
                        stop = Integer.MAX_VALUE;
                    } else {
                        start = parsePosition(str.substring(colonIndex + 1));   // chr1:1
                        stop = start;
                    }
                } else {
                    start = parsePosition(str.substring(colonIndex + 1, dashIndex));  // chr1:1-1
                    stop = parsePosition(str.substring(dashIndex + 1));
                }
            } catch(Exception e) {
                throw new UserException("Failed to parse Genome Location string: " + str, e);
            }
        }

        // is the contig valid?
        if (!isContigValid(contig))
            throw new UserException("Contig '" + contig + "' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?");

        if (stop == Integer.MAX_VALUE)
            // lookup the actually stop position!
            stop = getContigInfo(contig).getSequenceLength();

        GenomeLoc locus = new GenomeLoc(contig, getContigIndex(contig,true), start, stop);
        exceptionOnInvalidGenomeLoc(locus);
        return locus;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Parsing string representations
    //
    // --------------------------------------------------------------------------------------------------------------
    /**
     * Parses a number like 1,000,000 into a long.
     * @param pos
     */
    private int parsePosition(final String pos) {
        //String x = pos.replaceAll(",", "");  - this was replaced because it uses regexps
	//System.out.println("Parsing position: '" + pos + "'");
	if(pos.indexOf('-') != -1) {
	    throw new NumberFormatException("Position: '" + pos + "' can't contain '-'." );
	}

        if(pos.indexOf(',') != -1) {
            final StringBuilder buffer = new StringBuilder();
            for(int i = 0; i < pos.length(); i++) {
                final char c = pos.charAt(i);

                if(c == ',') {
                    continue;
                } else if(c < '0' || c > '9') {
                    throw new NumberFormatException("Position: '" + pos + "' contains invalid chars." );
                } else {
                    buffer.append(c);
                }
            }
            return Integer.parseInt(buffer.toString());
        } else {
            return Integer.parseInt(pos);
        }
    }


    /**
     * merge a list of genome locs that may be overlapping, returning the list of unique genomic locations
     *
     * @param raw the unchecked genome loc list
     * @param rule the merging rule we're using
     *
     * @return the list of merged locations
     */
    public List<GenomeLoc> mergeIntervalLocations(final List<GenomeLoc> raw, IntervalMergingRule rule) {
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

    /**
     * Determines whether the given contig is valid with respect to the sequence dictionary
     * already installed in the GenomeLoc.
     *
     * @return True if the contig is valid.  False otherwise.
     */
    private boolean isContigValid(String contig) {
        int contigIndex = contigInfo.getSequenceIndex(contig);
        return contigIndex >= 0 && contigIndex < contigInfo.size();
    }

    /**
     * Use this static constructor when the input data is under limited control (i.e. parsing user data).
     *
     * @param contig Contig to parse.
     * @param start  Starting point.
     * @param stop   Stop point.
     *
     * @return The genome location, or a MalformedGenomeLocException if unparseable.
     *
     * Validation: only checks that contig is valid
     * start/stop could be anything
     */
    public GenomeLoc parseGenomeLoc(final String contig, int start, int stop) {
        if (!isContigValid(contig))
            throw new MalformedGenomeLocException("Contig " + contig + " does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?");
        return new GenomeLoc(contig, getContigIndex(contig,true), start, stop);
    }


    /**
     * Read a file of genome locations to process.
     * regions specified by the location string.  The string is of the form:
     * Of the form: loc1;loc2;...
     * Where each locN can be:
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
     *
     * @param file_name interval file
     * @param allowEmptyIntervalList if false empty interval lists will return null
     * @return List<GenomeLoc> List of Genome Locs that have been parsed from file
     */
    public List<GenomeLoc> intervalFileToList(final String file_name, boolean allowEmptyIntervalList) {
        // try to open file
        File inputFile = new File(file_name);

        // check if file is empty
        if (inputFile.exists() && inputFile.length() < 1) {
            if (allowEmptyIntervalList)
                return new ArrayList<GenomeLoc>();
            else {
                Utils.warnUser("The interval file " + file_name + " is empty. The GATK will continue processing but you " +
                        "may want to fix (or exclude) this file.");
                return null;
            }
        }

        // case: BED file
        if (file_name.toUpperCase().endsWith(".BED")) {
            BedParser parser = new BedParser(this,inputFile);
            return parser.getLocations();
        }

        /**
         * IF not a BED file:
         * first try to read it as an interval file since that's well structured
         * we'll fail quickly if it's not a valid file.  Then try to parse it as
         * a location string file
         */
        try {
            IntervalList il = IntervalList.fromFile(inputFile);

            // iterate through the list of merged intervals and add then as GenomeLocs
            List<GenomeLoc> ret = new ArrayList<GenomeLoc>();
            for (Interval interval : il.getUniqueIntervals()) {
                ret.add(new GenomeLoc(interval.getSequence(), getContigIndex(interval.getSequence(),true), interval.getStart(), interval.getEnd()));
            }
            // always return null instead of empty list
            return ret.isEmpty() ? null : ret;

        }

        // if that didn't work, try parsing file as an old fashioned string file
        catch (Exception e) {
            try {
                List<GenomeLoc> ret = new ArrayList<GenomeLoc>();
                XReadLines reader = new XReadLines(new File(file_name));
                for(String line: reader) {
                        ret.add(parseGenomeInterval(line));
                }
                reader.close();

                // always return null instead of empty list
                return ret.isEmpty() ? null : ret;
            }
            catch (IOException e2) {
                throw new UserException.CouldNotReadInputFile(new File(file_name), e);
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
    private String getSequenceNameFromIndex(int contigIndex) {
        return contigInfo.getSequence(contigIndex).getSequenceName();
    }

    /**
     * create a genome loc, given the contig name, start, and stop
     *
     * @param contig the contig name
     * @param start  the starting position
     * @param stop   the stop position
     *
     * @return a new genome loc - but don't exception out if it is invalid
     */
    public GenomeLoc createPotentiallyInvalidGenomeLoc(String contig, final int start, final int stop) {
        return new GenomeLoc(contig, getContigIndex(contig,false), start, stop);
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
    public GenomeLoc createGenomeLoc(String contig, final int start, final int stop) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contig, getContigIndex(contig,true), start, stop));
    }

    /**
     * create a genome loc, given a read
     *
     * @param read
     *
     * @return
     */
    public GenomeLoc createGenomeLoc(final SAMRecord read) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(read.getReferenceName(), read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd()));
    }

    /**
     * create a new genome loc, given the contig name, and a single position
     *
     * @param contig the contig name
     * @param pos    the postion
     *
     * @return a genome loc representing a single base at the specified postion on the contig
     */
    public GenomeLoc createGenomeLoc(final String contig, final int pos) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contig, getContigIndex(contig,true), pos, pos));
    }

    /**
     * verify the specified genome loc is valid, if it's not, throw an exception
     * Will not verify the location against contig bounds.
     *
     *
     * Validation:
     * checks that start and stop are positive, start < stop, and the contig is valid
     * does not check that genomeLoc is actually on the contig, so start could be > end of contig
     *
     * @param toReturn the genome loc we're about to return
     *
     * @return the genome loc if it's valid, otherwise we throw an exception
     *
     */
    private GenomeLoc exceptionOnInvalidGenomeLoc(GenomeLoc toReturn) {
        if (toReturn.getStart() < 0) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the start position is less than 0");
        }
        if ((toReturn.getStop() != -1) && (toReturn.getStop() < 0)) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the stop position is less than 0");
        }
        if (toReturn.getContigIndex() < 0) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the contig index is less than 0");
        }
        if (toReturn.getContigIndex() >= contigInfo.getSequences().size()) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the contig index is greater then the stored sequence count");

        }
        return toReturn;

    }

    /**
     * Verify the locus against the bounds of the contig.
     *
     * performs boundary validation for genome loc INTERVALS:
     * start and stop are on contig and start <= stop
     * does NOT check that start and stop > 0, or that contig is valid
     * for that reason, this function should only be called AFTER exceptionOnInvalidGenomeLoc()
     * exceptionOnInvalidGenomeLoc isn't included in this function to save time
     *
     * @param locus Locus to verify.
     */
    private void exceptionOnInvalidGenomeLocBounds(GenomeLoc locus) {
        int contigSize = contigInfo.getSequence(locus.getContigIndex()).getSequenceLength();
        if(locus.getStart() > contigSize)
            throw new ReviewedStingException(String.format("GenomeLoc is invalid: locus start %d is after the end of contig %s",locus.getStart(),locus.getContig()));
        if(locus.getStop() > contigSize)
            throw new ReviewedStingException(String.format("GenomeLoc is invalid: locus stop %d is after the end of contig %s",locus.getStop(),locus.getContig()));
        if (locus.getStart() > locus.getStop()) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the start position is greater than the end position");
        }
    }

    /**
     * a method for validating genome locs as valid
     *
     * @param loc the location to validate
     *
     * @return true if the passed in GenomeLoc represents a valid location
     *
     * performs interval-style validation: contig is valid and atart and stop less than the end
     */
    public boolean validGenomeLoc(GenomeLoc loc) {
        // quick check before we get the contig size, is the contig number valid
        if ((loc.getContigIndex() < 0) ||                                       // the contig index has to be positive
            (loc.getContigIndex() >= contigInfo.getSequences().size()))         // the contig must be in the integer range of contigs)
            return false;

        int contigSize = contigInfo.getSequence(loc.getContigIndex()).getSequenceLength();
        if ((loc.getStart() < 0) ||                                             // start must be greater than 0
            ((loc.getStop() != -1) && (loc.getStop() < 0)) ||                   // the stop can be -1, but no other neg number
            (loc.getStart() > contigSize) ||                                    // the start must be before or equal to the contig end
            (loc.getStop() > contigSize))                                       // the stop must also be before or equal to the contig end
            return false;

        // we passed
        return true;
    }

    /**
     * validate a position or interval on the genome as valid
     *
     * @param contig the contig name
     * @param start  the start position
     * @param stop   the stop position
     *
     * @return true if it's valid, false otherwise
     *
     * performs interval-style validation: contig is valid and atart and stop less than the end
     */
    public boolean validGenomeLoc(String contig, int start, int stop) {
        return validGenomeLoc(new GenomeLoc(contig, getContigIndex(contig, false), start, stop));

    }

    /**
     * validate a position or interval on the genome as valid
     *
     * @param contigIndex the contig name
     * @param start       the start position
     * @param stop        the stop position
     *
     * @return true if it's valid, false otherwise
     *
     * performs interval-style validation: contig is valid and atart and stop less than the end
     */
    public boolean validGenomeLoc(int contigIndex, int start, int stop) {
        if (contigIndex < 0 || contigIndex >= contigInfo.size()) return false;
        return validGenomeLoc(new GenomeLoc(getSequenceNameFromIndex(contigIndex), contigIndex, start, stop));
    }

    /**
     * create a new genome loc from an existing loc, with a new start position
     * Note that this function will NOT explicitly check the ending offset, in case someone wants to
     * set the start of a new GenomeLoc pertaining to a read that goes off the end of the contig.
     *
     * @param loc   the old location
     * @param start a new start position
     *
     * @return the newly created genome loc
     */
    public GenomeLoc setStart(GenomeLoc loc, int start) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), start, loc.getStop()));
    }

    /**
     * create a new genome loc from an existing loc, with a new stop position
     * Note that this function will NOT explicitly check the ending offset, in case someone wants to
     * set the stop of a new GenomeLoc pertaining to a read that goes off the end of the contig.
     *
     * @param loc  the old location
     * @param stop a new stop position
     *
     * @return
     */
    public GenomeLoc setStop(GenomeLoc loc, int stop) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start, stop));
    }

    /**
     * return a new genome loc, with an incremented position
     *
     * @param loc the old location
     *
     * @return a new genome loc
     */
    public GenomeLoc incPos(GenomeLoc loc) {
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
    public GenomeLoc incPos(GenomeLoc loc, int by) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(loc.getContig(), loc.getContigIndex(), loc.start + by, loc.stop + by));
    }

    /**
     * Creates a GenomeLoc than spans the entire contig.
     * @param contigName Name of the contig.
     * @return A locus spanning the entire contig.
     */
    public GenomeLoc createOverEntireContig(String contigName) {
        SAMSequenceRecord contig = contigInfo.getSequence(contigName);
        if(contig == null)
            throw new ReviewedStingException("Unable to find contig named " + contigName);
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contigName,contig.getSequenceIndex(),1,contig.getSequenceLength()));
    }    
}
