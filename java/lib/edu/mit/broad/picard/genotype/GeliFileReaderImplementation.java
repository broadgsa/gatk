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
package edu.mit.broad.picard.genotype;


import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.io.StringReader;
import java.util.Arrays;

import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMTextHeaderCodec;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.BlockCompressedInputStream;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.StringLineReader;

/**
 * Internal class for reading GELI files.
 */
class GeliFileReaderImplementation extends GeliFileReader.ReaderImplementation {
    
    private boolean mIsSeekable = false;
    private BinaryCodec mStream = null;
    private final BlockCompressedInputStream mCompressedInputStream;
    private SAMFileHeader mFileHeader = null;
    private long mFirstRecordPointer = 0;
    private CloseableIterator<GenotypeLikelihoods> mCurrentIterator = null;


    GeliFileReaderImplementation(final InputStream stream)
        throws IOException {
        mIsSeekable = false;
        mCompressedInputStream = new BlockCompressedInputStream(stream);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        readHeader(null);
    }

    GeliFileReaderImplementation(final File file)
        throws IOException {
        mIsSeekable = true;
        mCompressedInputStream = new BlockCompressedInputStream(file);
        mStream = new BinaryCodec(new DataInputStream(mCompressedInputStream));
        readHeader(file);
        mFirstRecordPointer = mCompressedInputStream.getFilePointer();
    }

    void close() {
        if (mStream != null) {
            mStream.close();
        }
        mStream = null;
        mFileHeader = null;
    }

    SAMFileHeader getFileHeader() {
        return mFileHeader;
    }

    CloseableIterator<GenotypeLikelihoods> getIterator() {
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
        mCurrentIterator = new GELIFileIterator();
        return mCurrentIterator;
    }

    private void readHeader(final File file)
        throws IOException {

        final byte[] buffer = new byte[4];
        mStream.readBytes(buffer);
        if (!Arrays.equals(buffer, GeliFileConstants.GELI_MAGIC)) {
            throw new IOException("Invalid GELI file header");
        }

        final int headerTextLength = mStream.readInt();
        final String textHeader = mStream.readString(headerTextLength);
        mFileHeader = new SAMTextHeaderCodec().decode(new StringLineReader(textHeader),
                file);

        final int sequenceCount = mStream.readInt();
        if (sequenceCount != mFileHeader.getSequences().size()) {
            throw new GeliException("Number of sequences in text header (" + mFileHeader.getSequences().size() +
            ") != number of sequences in binary header (" + sequenceCount + ") for file " + file);
        }
        for (int i = 0; i < sequenceCount; i++) {
            readSequenceRecord(file);
//            final SAMSequenceRecord sequenceRecord = mFileHeader.getSequence(i);
//            if (!sequenceRecord.getSequenceName().equals(binarySequenceRecord.getSequenceName())) {
//                throw new GELIException("For sequence " + i + ", text and binary have different names in file " +
//                file);
//            }
//            if (sequenceRecord.getSequenceLength() != binarySequenceRecord.getSequenceLength()) {
//                throw new GELIException("For sequence " + i + ", text and binary have different lengths in file " +
//                file);
//            }
        }
    }

    private SAMSequenceRecord readSequenceRecord(final File file) {
        final int nameLength = mStream.readInt();
        if (nameLength <= 1) {
            throw new GeliException("Invalid BAM file header: missing sequence name in file " + file);
        }
        final String sequenceName = mStream.readString(nameLength - 1);
        // Skip the null terminator
        mStream.readByte();
        final int sequenceLength = mStream.readInt();
        final SAMSequenceRecord record = new SAMSequenceRecord(sequenceName);
        record.setSequenceLength(sequenceLength);
        return record;
    }

    private class GELIFileIterator
        implements CloseableIterator<GenotypeLikelihoods> {

        private GenotypeLikelihoods mNextRecord = null;
        private final GenotypeLikelihoodsCodec likelihoodsCodec = new GenotypeLikelihoodsCodec();


        GELIFileIterator() {
            this(true);
        }

        GELIFileIterator(final boolean advance) {
            likelihoodsCodec.setInputStream(mStream.getInputStream());
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

        public GenotypeLikelihoods next() {
            final GenotypeLikelihoods result = mNextRecord;
            advance();
            return result;
        }

        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }

        void advance() {
            try {
                mNextRecord = getNextRecord();
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
        }

        GenotypeLikelihoods getNextRecord()
            throws IOException {
            return likelihoodsCodec.decode();
        }
    }
}
