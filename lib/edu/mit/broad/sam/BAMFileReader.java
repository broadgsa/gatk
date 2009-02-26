/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.sam;


import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.BlockCompressedInputStream;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.StringLineReader;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Internal class for reading and querying BAM files.
 */
class BAMFileReader
    extends SAMFileReader.ReaderImplementation {
    private boolean mIsSeekable = false;
    private BinaryCodec mStream = null;
    private final BlockCompressedInputStream mCompressedInputStream;
    private SAMFileHeader mFileHeader = null;
    private BAMFileIndex mFileIndex = null;
    private long mFirstRecordPointer = 0;
    private CloseableIterator<SAMRecord> mCurrentIterator = null;
    private final boolean eagerDecode;


    BAMFileReader(final InputStream stream, final boolean eagerDecode)
        throws IOException {
        mIsSeekable = false;
        mCompressedInputStream = new BlockCompressedInputStream(stream);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        readHeader(null);
    }

    BAMFileReader(final File file, final boolean eagerDecode)
        throws IOException {
        mIsSeekable = true;
        mCompressedInputStream = new BlockCompressedInputStream(file);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        this.eagerDecode = eagerDecode;
        readHeader(file);
        mFirstRecordPointer = mCompressedInputStream.getFilePointer();
    }

    void close() {
        if (mStream != null) {
            mStream.close();
        }
        mStream = null;
        mFileHeader = null;
        mFileIndex = null;
    }

    BAMFileIndex getFileIndex() {
        return mFileIndex;
    }

    void setFileIndex(final BAMFileIndex fileIndex) {
        mFileIndex = fileIndex;
    }

    SAMFileHeader getFileHeader() {
        return mFileHeader;
    }

    /**
     * Currently this is ignored for BAM reading.  Always do strict validation.
     */
    void setValidationStringency(final SAMFileReader.ValidationStringency validationStringency) {
    }

    CloseableIterator<SAMRecord> getIterator() {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (mIsSeekable) {
            try {
                mCompressedInputStream.seek(mFirstRecordPointer);
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
        }
        mCurrentIterator = new BAMFileIterator();
        return mCurrentIterator;
    }

    CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        if (mStream == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        if (!mIsSeekable) {
            throw new UnsupportedOperationException("Cannot query stream-based BAM file");
        }
        if (mFileIndex == null) {
            throw new IllegalStateException("No BAM file index is available");
        }
        mCurrentIterator = new BAMFileIndexIterator(sequence, start, end, contained);
        return mCurrentIterator;
    }

    private void readHeader(final File file)
        throws IOException {

        final byte[] buffer = new byte[4];
        mStream.readBytes(buffer);
        if (!Arrays.equals(buffer, BAMFileConstants.BAM_MAGIC)) {
            throw new IOException("Invalid BAM file header");
        }

        final int headerTextLength = mStream.readInt();
        final String textHeader = mStream.readString(headerTextLength);
        mFileHeader = new SAMTextHeaderCodec().decode(new StringLineReader(textHeader),
                file);

        final int sequenceCount = mStream.readInt();
        if (mFileHeader.getSequences().size() > 0) {
            // It is allowed to have binary sequences but no text sequences, so only validate if both are present
            if (sequenceCount != mFileHeader.getSequences().size()) {
                throw new SAMFormatException("Number of sequences in text header (" + mFileHeader.getSequences().size() +
                        ") != number of sequences in binary header (" + sequenceCount + ") for file " + file);
            }
            for (int i = 0; i < sequenceCount; i++) {
                final SAMSequenceRecord binarySequenceRecord = readSequenceRecord(file);
                final SAMSequenceRecord sequenceRecord = mFileHeader.getSequence(i);
                if (!sequenceRecord.getSequenceName().equals(binarySequenceRecord.getSequenceName())) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different names in file " +
                            file);
                }
                if (sequenceRecord.getSequenceLength() != binarySequenceRecord.getSequenceLength()) {
                    throw new SAMFormatException("For sequence " + i + ", text and binary have different lengths in file " +
                            file);
                }
            }
        } else {
            // If only binary sequences are present, copy them into mFileHeader
            final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>(sequenceCount);
            for (int i = 0; i < sequenceCount; i++) {
                sequences.add(readSequenceRecord(file));
            }
            mFileHeader.setSequences(sequences);
        }
    }

    private SAMSequenceRecord readSequenceRecord(final File file) {
        final int nameLength = mStream.readInt();
        if (nameLength <= 1) {
            throw new SAMFormatException("Invalid BAM file header: missing sequence name in file " + file);
        }
        final String sequenceName = mStream.readString(nameLength - 1);
        // Skip the null terminator
        mStream.readByte();
        final int sequenceLength = mStream.readInt();
        final SAMSequenceRecord record = new SAMSequenceRecord(sequenceName);
        record.setSequenceLength(sequenceLength);
        return record;
    }

    private class BAMFileIterator
        implements CloseableIterator<SAMRecord> {

        private SAMRecord mNextRecord = null;
        private final BAMRecordCodec bamRecordCodec = new BAMRecordCodec(getFileHeader());


        BAMFileIterator() {
            this(true);
        }

        BAMFileIterator(final boolean advance) {
            this.bamRecordCodec.setInputStream(BAMFileReader.this.mStream.getInputStream());

            if (advance) {
                advance();
            }
        }

        public void close() {
            if (this != mCurrentIterator) {
                throw new IllegalStateException("Attempt to close non-current iterator");
            }
            mCurrentIterator = null;
        }

        public boolean hasNext() {
            return (mNextRecord != null);
        }

        public SAMRecord next() {
            final SAMRecord result = mNextRecord;
            advance();
            return result;
        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }

        void advance() {
            try {
                mNextRecord = getNextRecord();
                if (eagerDecode && mNextRecord != null) {
                    mNextRecord.eagerDecode();
                }
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
        }

        SAMRecord getNextRecord()
            throws IOException {
            return bamRecordCodec.decode();
        }
    }

    private class BAMFileIndexIterator
        extends BAMFileIterator {

        private long[] mFilePointers = null;
        private int mFilePointerIndex = 0;
        private long mFilePointerLimit = -1;
        private int mReferenceIndex = -1;
        private int mRegionStart = 0;
        private int mRegionEnd = 0;
        private boolean mReturnContained = false;


        BAMFileIndexIterator(final String sequence, final int start, final int end, final boolean contained) {
            super(false);  // delay advance() until after construction
            final SAMFileHeader fileHeader = getFileHeader();
            mReferenceIndex = fileHeader.getSequenceIndex(sequence);
            if (mReferenceIndex != -1) {
                final BAMFileIndex fileIndex = getFileIndex();
                mFilePointers = fileIndex.getSearchBins(mReferenceIndex, start, end);
            }
            mRegionStart = start;
            mRegionEnd = (end <= 0) ? Integer.MAX_VALUE : end;
            mReturnContained = contained;
            advance();
        }

        SAMRecord getNextRecord()
            throws IOException {
            while (true) {
                // Advance to next file block if necessary
                while (mCompressedInputStream.getFilePointer() >= mFilePointerLimit) {
                    if (mFilePointers == null ||
                        mFilePointerIndex >= mFilePointers.length) {
                        return null;
                    }
                    final long startOffset = mFilePointers[mFilePointerIndex++];
                    final long endOffset = mFilePointers[mFilePointerIndex++];
                    mCompressedInputStream.seek(startOffset);
                    mFilePointerLimit = endOffset;
                }
                // Pull next record from stream
                final SAMRecord record = super.getNextRecord();
                if (record == null) {
                    return null;
                }
                // If beyond the end of this reference sequence, end iteration
                final int referenceIndex = record.getReferenceIndex();
                if (referenceIndex != mReferenceIndex) {
                    if (referenceIndex < 0 ||
                        referenceIndex > mReferenceIndex) {
                        mFilePointers = null;
                        return null;
                    }
                    // If before this reference sequence, continue
                    continue;
                }
                if (mRegionStart == 0 && mRegionEnd == Integer.MAX_VALUE) {
                    // Quick exit to avoid expensive alignment end calculation
                    return record;
                }
                final int alignmentStart = record.getAlignmentStart();
                final int alignmentEnd = record.getAlignmentEnd();
                if (alignmentStart > mRegionEnd) {
                    // If scanned beyond target region, end iteration
                    mFilePointers = null;
                    return null;
                }
                // Filter for overlap with region
                if (mReturnContained) {
                    if (alignmentStart >= mRegionStart && alignmentEnd <= mRegionEnd) {
                        return record;
                    }
                } else {
                    if (alignmentEnd >= mRegionStart && alignmentStart <= mRegionEnd) {
                        return record;
                    }
                }
            }
        }
    }
}
