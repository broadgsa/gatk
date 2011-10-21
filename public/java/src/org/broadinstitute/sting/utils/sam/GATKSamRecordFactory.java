/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordFactory;
import net.sf.samtools.BAMRecord;
import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Factory interface implementation used to create GATKSamRecords
 * from SAMFileReaders with SAM-JDK
 *
 * @author Mark DePristo
 */
public class GATKSamRecordFactory implements SAMRecordFactory {

    /** Create a new SAMRecord to be filled in */
    public SAMRecord createSAMRecord(SAMFileHeader header) {
        throw new UserException.BadInput("The GATK now longer supports input SAM files");
    }

    /** Create a new BAM Record. */
    public BAMRecord createBAMRecord(final SAMFileHeader header,
                                     final int referenceSequenceIndex,
                                     final int alignmentStart,
                                     final short readNameLength,
                                     final short mappingQuality,
                                     final int indexingBin,
                                     final int cigarLen,
                                     final int flags,
                                     final int readLen,
                                     final int mateReferenceSequenceIndex,
                                     final int mateAlignmentStart,
                                     final int insertSize,
                                     final byte[] variableLengthBlock) {
        return new GATKSAMRecord(header,
                referenceSequenceIndex,
                alignmentStart,
                readNameLength,
                mappingQuality,
                indexingBin,
                cigarLen,
                flags,
                readLen,
                mateReferenceSequenceIndex,
                mateAlignmentStart,
                insertSize,
                variableLengthBlock);
    }
}
