/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.filters;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Filter out malformed reads.
 *
 * @author mhanna
 * @version 0.1
 */
public class MalformedReadFilter extends ReadFilter {
    private SAMFileHeader header;

    @Argument(fullName = "filter_mismatching_base_and_quals", shortName = "filterMBQ", doc = "if a read has mismatching number of bases and base qualities, filter out the read instead of blowing up.", required = false)
    boolean filterMismatchingBaseAndQuals = false;

    @Override
    public void initialize(GenomeAnalysisEngine engine) {
        this.header = engine.getSAMFileHeader();
    }

    public boolean filterOut(SAMRecord read) {
        // slowly changing the behavior to blow up first and filtering out if a parameter is explicitly provided
        if (!checkMismatchingBasesAndQuals(read)) {
            if (!filterMismatchingBaseAndQuals)
                throw new UserException.MalformedBAM(read, "BAM file has a read with mismatching number of bases and base qualities. Offender: " + read.getReadName() +"  [" + read.getReadLength() + " bases] [" +read.getBaseQualities().length +"] quals");
            else
                return true;
        }

        return  !checkInvalidAlignmentStart(read) ||
                !checkInvalidAlignmentEnd(read) ||
                !checkAlignmentDisagreesWithHeader(this.header,read) ||
                !checkCigarDisagreesWithAlignment(read);
    }

    /**
     * Check for the case in which the alignment start is inconsistent with the read unmapped flag.
     * @param read The read to validate.
     * @return true if read start is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentStart( SAMRecord read ) {
        // read is not flagged as 'unmapped', but alignment start is NO_ALIGNMENT_START
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START )
            return false;
        // Read is not flagged as 'unmapped', but alignment start is -1
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == -1 )
            return false;
        return true;
    }

    /**
     * Check for invalid end of alignments.
     * @param read The read to validate.
     * @return true if read end is valid, false otherwise.
     */
    private static boolean checkInvalidAlignmentEnd( SAMRecord read ) {
        // Alignment aligns to negative number of bases in the reference.
        if( !read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1 && (read.getAlignmentEnd()-read.getAlignmentStart()+1)<0 )
            return false;
        return true;
    }

    /**
     * Check to ensure that the alignment makes sense based on the contents of the header.
     * @param header The SAM file header.
     * @param read The read to verify.
     * @return true if alignment agrees with header, false othrewise.
     */
    private static boolean checkAlignmentDisagreesWithHeader( SAMFileHeader header, SAMRecord read ) {
        // Read is aligned to nonexistent contig
        if( read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START )
            return false;
        SAMSequenceRecord contigHeader = header.getSequence( read.getReferenceIndex() );
        // Read is aligned to a point after the end of the contig
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() > contigHeader.getSequenceLength() )
            return false;
        return true;
    }

    /**
     * Check for inconsistencies between the cigar string and the
     * @param read The read to validate.
     * @return true if cigar agrees with alignment, false otherwise.
     */
    private static boolean checkCigarDisagreesWithAlignment(SAMRecord read) {
        // Read has a valid alignment start, but the CIGAR string is empty
        if( !read.getReadUnmappedFlag() &&
            read.getAlignmentStart() != -1 &&
            read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START &&
            read.getAlignmentBlocks().size() < 0 )
            return false;
        return true;
    }

    /**
     * Check if the read has the same number of bases and base qualities
     * @param read the read to validate
     * @return true if they have the same number. False otherwise.
     */
    private static boolean checkMismatchingBasesAndQuals(SAMRecord read) {
        return (read.getReadLength() == read.getBaseQualities().length);
    }
}
