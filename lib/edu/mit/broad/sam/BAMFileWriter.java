/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.BlockCompressedOutputStream;

import java.io.DataOutputStream;
import java.io.File;

/**
 * Concrete implementation of SAMFileWriter for writing gzipped BAM files.
 */
class BAMFileWriter extends SAMFileWriterImpl {

    private final BinaryCodec outputBinaryCodec;
    private BAMRecordCodec bamRecordCodec = null;

    public BAMFileWriter(final File path) {
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(new BlockCompressedOutputStream(path)));
        outputBinaryCodec.setOutputFileName(path.toString());
    }

    private void prepareToWriteAlignments() {
        if (bamRecordCodec == null) {
            bamRecordCodec = new BAMRecordCodec(getHeader());
            bamRecordCodec.setOutputStream(outputBinaryCodec.getOutputStream());
        }
    }

    protected void writeAlignment(final SAMRecord alignment) {
        prepareToWriteAlignments();
        bamRecordCodec.encode(alignment);
    }

    protected void writeHeader(final String textHeader) {
        outputBinaryCodec.writeBytes(BAMFileConstants.BAM_MAGIC);

        // calculate and write the length of the SAM file header text and the header text
        outputBinaryCodec.writeString(textHeader, true, false);

        // write the sequences binarily.  This is redundant with the text header
        outputBinaryCodec.writeInt(getHeader().getSequences().size());
        for (final SAMSequenceRecord sequenceRecord: getHeader().getSequences()) {
            outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }
    }

    protected void finish() {
        outputBinaryCodec.close();
    }

    protected String getFilename() {
        return outputBinaryCodec.getOutputFileName();
    }
}
