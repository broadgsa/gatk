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

import java.util.List;

import com.google.java.contract.*;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Factory class for creating GenomeLocs
 */
@Invariant({
        "logger != null",
        "contigInfo != null"})
public class GenomeLocParser {
    private static Logger logger = Logger.getLogger(GenomeLocParser.class);

    // --------------------------------------------------------------------------------------------------------------
    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    // --------------------------------------------------------------------------------------------------------------
    private final MasterSequenceDictionary contigInfo;

    /**
     * A wrapper class that provides efficient last used caching for the global
     * SAMSequenceDictionary underlying all of the GATK engine capabilities
     */
    // todo -- enable when CoFoJa developers identify the problem (likely thread unsafe invariants)
//    @Invariant({
//            "dict != null",
//            "dict.size() > 0",
//            "lastSSR == null || dict.getSequence(lastContig).getSequenceIndex() == lastIndex",
//            "lastSSR == null || dict.getSequence(lastContig).getSequenceName() == lastContig",
//            "lastSSR == null || dict.getSequence(lastContig) == lastSSR"})
    private final class MasterSequenceDictionary {
        final private SAMSequenceDictionary dict;

        // cache
        SAMSequenceRecord lastSSR = null;
        String lastContig = "";
        int lastIndex = -1;

        @Requires({"dict != null", "dict.size() > 0"})
        public MasterSequenceDictionary(SAMSequenceDictionary dict) {
            this.dict = dict;
        }

        @Ensures("result > 0")
        public final int getNSequences() {
            return dict.size();
        }

        @Requires("contig != null")
        public synchronized boolean hasContig(final String contig) {
            return lastContig == contig || dict.getSequence(contig) != null;
        }

        @Requires("index >= 0")
        public synchronized boolean hasContig(final int index) {
            return lastIndex == index|| dict.getSequence(index) != null;
        }

        @Requires("contig != null")
        @Ensures("result != null")
        public synchronized final SAMSequenceRecord getSequence(final String contig) {
            if ( isCached(contig) )
                return lastSSR;
            else
                return updateCache(contig, -1);
        }

        @Requires("index >= 0")
        @Ensures("result != null")
        public synchronized final SAMSequenceRecord getSequence(final int index) {
            if ( isCached(index) )
                return lastSSR;
            else
                return updateCache(null, index);

        }

        @Requires("contig != null")
        @Ensures("result >= 0")
        public synchronized final int getSequenceIndex(final String contig) {
            if ( ! isCached(contig) ) {
                updateCache(contig, -1);
            }

            return lastIndex;
        }

        @Requires({"contig != null", "lastContig != null"})
        private synchronized boolean isCached(final String contig) {
            return lastContig.equals(contig);
        }

        @Requires({"lastIndex != -1", "index >= 0"})
        private synchronized boolean isCached(final int index) {
            return lastIndex == index;
        }

        /**
         * The key algorithm.  Given a new record, update the last used record, contig
         * name, and index.
         *
         * @param contig
         * @param index
         * @return
         */
        @Requires("contig != null || index >= 0")
        @Ensures("result != null")
        private synchronized SAMSequenceRecord updateCache(final String contig, int index ) {
            SAMSequenceRecord rec = contig == null ? dict.getSequence(index) : dict.getSequence(contig);
            if ( rec == null ) {
                throw new ReviewedStingException("BUG: requested unknown contig=" + contig + " index=" + index);
            } else {
                lastSSR = rec;
                lastContig = rec.getSequenceName();
                lastIndex = rec.getSequenceIndex();
                return rec;
            }
        }


    }

    /**
     * set our internal reference contig order
     * @param refFile the reference file
     */
    @Requires("refFile != null")
    public GenomeLocParser(final ReferenceSequenceFile refFile) {
        this(refFile.getSequenceDictionary());
    }

    public GenomeLocParser(SAMSequenceDictionary seqDict) {
        if (seqDict == null) { // we couldn't load the reference dictionary
            //logger.info("Failed to load reference dictionary, falling back to lexicographic order for contigs");
            throw new UserException.CommandLineException("Failed to load reference dictionary");
        }

        contigInfo = new MasterSequenceDictionary(seqDict);
        logger.debug(String.format("Prepared reference sequence contig dictionary"));
        for (SAMSequenceRecord contig : seqDict.getSequences()) {
            logger.debug(String.format(" %s (%d bp)", contig.getSequenceName(), contig.getSequenceLength()));
        }
    }

    /**
     * get the contig's SAMSequenceRecord
     *
     * @param contig the string name of the contig
     *
     * @return the sam sequence record
     */
    @Requires({"contig != null", "contigIsInDictionary(contig)"})
    @Ensures("result != null")
    public SAMSequenceRecord getContigInfo(final String contig) {
        return contigInfo.getSequence(contig);
    }

    /**
     * Determines whether the given contig is valid with respect to the sequence dictionary
     * already installed in the GenomeLoc.
     *
     * @return True if the contig is valid.  False otherwise.
     */
    @Requires({"contig != null"})
    public boolean contigIsInDictionary(String contig) {
        return contigInfo.hasContig(contig);
    }

    public boolean indexIsInDictionary(final int index) {
        return contigInfo.hasContig(index);
    }

    /**
     * Returns the contig index of a specified string version of the contig
     *
     * @param contig the contig string
     *
     * @return the contig index, -1 if not found
     */
    @Requires("contig != null")
    @Ensures("result >= 0")
    public int getContigIndex(final String contig) {
        if (! contigInfo.hasContig(contig) )
            throw new UserException.MalformedGenomeLoc(String.format("Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig));
        return contigInfo.getSequenceIndex(contig);
    }

    @Requires("contig != null")
    protected int getContigIndexWithoutException(final String contig) {
        if (! contigInfo.hasContig(contig))
            return -1;
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
    @Requires("str != null")
    @Ensures("result != null")
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
    @Requires("str != null")
    @Ensures("result != null")
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
        if (!contigIsInDictionary(contig))
            throw new UserException("Contig '" + contig + "' does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?");

        if (stop == Integer.MAX_VALUE)
            // lookup the actually stop position!
            stop = getContigInfo(contig).getSequenceLength();

        GenomeLoc locus = new GenomeLoc(contig, getContigIndex(contig), start, stop);
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
    @Requires("pos != null")
    @Ensures("result >= 0")
    private int parsePosition(final String pos) {
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
        if (!contigIsInDictionary(contig))
            throw new MalformedGenomeLocException("Contig " + contig + " does not match any contig in the GATK sequence dictionary derived from the reference; are you sure you are using the correct reference fasta file?");
        return new GenomeLoc(contig, getContigIndex(contig), start, stop);
    }

    /**
     * get the sequence name from a sequence index
     *
     * @param contigIndex get the contig index
     *
     * @return the string that represents that contig name
     */
    @Requires("indexIsInDictionary(contigIndex)")
    @Ensures("result != null")
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
     * @return a new genome loc
     */
    public GenomeLoc createGenomeLoc(String contig, final int start, final int stop) {
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contig, getContigIndex(contig), start, stop));
    }

    /**
     * create a genome loc, given a read. If the read is unmapped, *and* yet the read has a contig and start position,
     * then a GenomeLoc is returned for contig:start-start, otherwise and UNMAPPED GenomeLoc is returned.
     *
     * @param read
     *
     * @return
     */
    public GenomeLoc createGenomeLoc(final SAMRecord read) {
        GenomeLoc loc;

        if ( read.getReadUnmappedFlag() && read.getReferenceIndex() == -1 )
            // read is unmapped and not placed anywhere on the genome
            loc = GenomeLoc.UNMAPPED;
        else {
            int end = read.getReadUnmappedFlag() ? read.getAlignmentStart() : read.getAlignmentEnd();
            loc = new GenomeLoc(read.getReferenceName(), read.getReferenceIndex(), read.getAlignmentStart(), end);
        }

        return exceptionOnInvalidGenomeLoc(loc);
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
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contig, getContigIndex(contig), pos, pos));
    }

    /**
     * Creates a new GenomeLoc without performing any validation on its contig or bounds.
     * FOR UNIT TESTING PURPOSES ONLY!
     *
     * @param contig the contig name
     * @param start  start position of the interval
     * @param stop   stop position of the interval
     * @return a new GenomeLoc representing the specified location
     */
    public GenomeLoc createGenomeLocWithoutValidation( String contig, int start, int stop ) {
        return new GenomeLoc(contig, getContigIndexWithoutException(contig), start, stop);
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
        if ( GenomeLoc.isUnmapped(toReturn) ) {
            return toReturn;
        }

        if (toReturn.getStart() < 0) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the start position is less than 0 " +
                                             "in interval: " + toReturn);
        }
        if ((toReturn.getStop() != -1) && (toReturn.getStop() < 0)) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the stop position is less than 0 " +
                                             "in interval: " + toReturn);
        }
        if (toReturn.getContigIndex() < 0) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the contig index is less than 0 " +
                                             "in interval: " + toReturn);
        }
        if (toReturn.getContigIndex() >= contigInfo.getNSequences()) {
            throw new ReviewedStingException("Parameters to GenomeLocParser are incorrect: the contig index is greater than " +
                                             "the stored sequence count in interval: " + toReturn);

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
        if ( GenomeLoc.isUnmapped(locus) ) {
            return;
        }

        int contigSize = contigInfo.getSequence(locus.getContigIndex()).getSequenceLength();
        if(locus.getStart() > contigSize)
            throw new UserException.MalformedGenomeLoc("GenomeLoc is invalid: locus start is after the end of contig",locus);
        if(locus.getStop() > contigSize)
            throw new UserException.MalformedGenomeLoc("GenomeLoc is invalid: locus stop is after the end of contig",locus);
        if (locus.getStart() > locus.getStop()) {
            throw new UserException.MalformedGenomeLoc("Parameters to GenomeLocParser are incorrect: the start position is " +
                                                       "greater than the end position", locus);
        }
    }

    /**
     * Validates each of the GenomeLocs in the list passed as an argument: each GenomeLoc must refer to a
     * valid contig, and its bounds must be within the bounds of its contig. Throws a MalformedGenomeLoc
     * exception if an invalid GenomeLoc is found.
     *
     * @param genomeLocList The list of GenomeLocs to validate
     */
    public void validateGenomeLocList( List<GenomeLoc> genomeLocList ) {
        if ( genomeLocList == null ) {
            return;
        }

        for ( GenomeLoc loc : genomeLocList ) {
            try {
                exceptionOnInvalidGenomeLoc(loc);
                exceptionOnInvalidGenomeLocBounds(loc);
            }
            catch ( UserException e ) {
                throw e;
            }
            catch ( ReviewedStingException e ) {
                throw new UserException.MalformedGenomeLoc(e.getMessage());
            }
        }
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
        if ( ! contigInfo.hasContig(contig) )
            return false;
        else
            return validGenomeLoc(getContigIndex(contig), start, stop);
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
        // quick check before we get the contig size, is the contig number valid
        if ((contigIndex < 0) ||                                       // the contig index has to be positive
            (contigIndex >= contigInfo.getNSequences()))               // the contig must be in the integer range of contigs)
            return false;

        int contigSize = contigInfo.getSequence(contigIndex).getSequenceLength();
        if ((start < 0) ||                                             // start must be greater than 0
            (stop < 0 ) ||                                             // the stop must be greater than 0
            (start > contigSize) ||                                    // the start must be before or equal to the contig end
            (stop > contigSize))                                       // the stop must also be before or equal to the contig end
            return false;

        // we passed
        return true;
    }

    /**
     * create a new genome loc from an existing loc, with a new start position
     * Note that this function will NOT explicitly check the ending offset, in case someone wants to
     * set the start of a new GenomeLoc pertaining to a read that goes off the end of the contig.
     *
     * @param loc   the old location
     * @param start a new start position
     * TODO -- move to me GenomeLoc class itself
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
     * TODO -- move to me GenomeLoc class itself
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
     * TODO -- move to me GenomeLoc class itself
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
     * TODO -- move to me GenomeLoc class itself
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
    @Requires("contigName != null")
    @Ensures("result != null")
    public GenomeLoc createOverEntireContig(String contigName) {
        SAMSequenceRecord contig = contigInfo.getSequence(contigName);
        return exceptionOnInvalidGenomeLoc(new GenomeLoc(contigName,contig.getSequenceIndex(),1,contig.getSequenceLength()));
    }    
}
