/*
 * Copyright (c) 2009 The Broad Institute
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

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Validates reads against a specific set of criteria.  If it finds a
 * read that fails to meet the given criteria, it will throw an exception.
 * The caller can decide whether to ignore the error, hide the read
 * from the user, or blow up in a spectacular ball of fire.
 *
 * @author hanna
 * @version 0.1
 */
public class SAMReadValidator {
    /**
     * Validate the sam read against a list of criteria that are known to cause failures in the GATK.
     * Throw an exception if the read fails.
     * @param read the read to validate.  Must not be null.
     */
    public static void validate( SAMFileHeader header, SAMRecord read ) throws SAMReadValidationException {
        checkInvalidAlignmentStart(read);
        checkInvalidAlignmentEnd(read);
        checkAlignmentDisagreesWithHeader(header,read);
        checkCigarDisagreesWithAlignment(read);
    }

    /**
     * Check for the case in which the alignment start is inconsistent with the read unmapped flag.
     * @param read The read to validate.
     */
    private static void checkInvalidAlignmentStart( SAMRecord read ) {
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START )
            throw new SAMReadValidationException("read is not flagged as 'unmapped', but alignment start is NO_ALIGNMENT_START");
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() == -1 )
            throw new SAMReadValidationException("Read is not flagged as 'unmapped', but alignment start is -1");
    }

    /**
     * Check for invalid end of alignments.
     * @param read The read to validate.
     */
    private static void checkInvalidAlignmentEnd( SAMRecord read ) {
        if( !read.getReadUnmappedFlag() && read.getAlignmentEnd() != -1 && read.getAlignmentEnd() < read.getAlignmentStart() )
            throw new SAMReadValidationException("Alignment ends prior to its beginning");
    }

    private static void checkAlignmentDisagreesWithHeader( SAMFileHeader header, SAMRecord read ) {
        SAMSequenceRecord contigHeader = header.getSequence( read.getReferenceIndex() );
        if( !read.getReadUnmappedFlag() && read.getAlignmentStart() > contigHeader.getSequenceLength() )
            throw new SAMReadValidationException("Read is aligned to a point after the end of the contig");
    }

    /**
     * Check for inconsistencies between the cigar string and the 
     * @param read The read to validate.
     */
    private static void checkCigarDisagreesWithAlignment( SAMRecord read ) {
        if( !read.getReadUnmappedFlag() &&
            read.getAlignmentStart() != -1 &&
            read.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START &&
            read.getAlignmentBlocks().size() == 0 )
            throw new SAMReadValidationException("Read has a valid alignment start, but the CIGAR string is empty");
    }
}

