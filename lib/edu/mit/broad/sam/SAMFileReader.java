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


import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.RuntimeIOException;
import edu.mit.broad.sam.util.BlockCompressedInputStream;

import java.io.*;


/**
 * Class for reading and querying SAM/BAM files.
 */
public class SAMFileReader implements Iterable<SAMRecord>
{
    private boolean mIsBinary = false;
    private BAMFileIndex mFileIndex = null;
    private ReaderImplementation mReader = null;

    public enum ValidationStringency {
        STRICT,  // Do the right thing, throw an exception if something looks wrong
        LENIENT, // Emit warnings but keep going if possible
        SILENT;  // Like LENIENT, only don't emit warning messages

        public static ValidationStringency DEFAULT_STRINGENCY = STRICT;
    }

    /**
     * Internal interface for SAM/BAM file reader implementations.
     * Implemented as an abstract class to enforce better access control.
     */
    static abstract class ReaderImplementation {
        abstract SAMFileHeader getFileHeader();
        abstract CloseableIterator<SAMRecord> getIterator();
        abstract CloseableIterator<SAMRecord> query(String sequence, int start, int end, boolean contained);
        abstract void close();
        // If true, emit warnings about format errors rather than throwing exceptions;
        abstract void setValidationStringency(final ValidationStringency validationStringency);
    }


    public SAMFileReader(final InputStream stream) {
        this(stream, false);
    }

    public SAMFileReader(final File file) {
        this(file, null, false);
    }

    public SAMFileReader(final File file, final File indexFile) {
        this(file, indexFile, false);
    }

    /**
     * Read a SAM or BAM file
     * @param stream input SAM or BAM
     * @param eagerDecode if true, decode SAM record entirely when reading it
     */
    public SAMFileReader(final InputStream stream, final boolean eagerDecode) {
        init(stream, eagerDecode);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file if present
     * @param file where to read from
     * @param eagerDecode if true, decode SAM record entirely when reading it
     */
    public SAMFileReader(final File file, final boolean eagerDecode) {
        init(file, null, eagerDecode);
    }

    /**
     * Read a SAM or BAM file, possibly with an index file
     * @param file where to read from
     * @param indexFile location of index file, or null in order to use the default index file (if present)
     * @param eagerDecode eagerDecode if true, decode SAM record entirely when reading it
     */
    public SAMFileReader(final File file, final File indexFile, final boolean eagerDecode){
        init(file, indexFile, eagerDecode);
    }

    public void close() {
        if (mReader != null) {
            mReader.close();
        }
        if (mFileIndex != null) {
            mFileIndex.close();
        }
        mReader = null;
        mFileIndex = null;
    }

    public boolean isBinary() {
        return mIsBinary;
    }

    public boolean hasIndex() {
        return (mFileIndex != null);
    }

    public SAMFileHeader getFileHeader() {
        return mReader.getFileHeader();
    }

    public void setValidationStringency(final ValidationStringency validationStringency) {
        mReader.setValidationStringency(validationStringency);
    }

    public CloseableIterator<SAMRecord> iterator() {
        return mReader.getIterator();
    }

    public CloseableIterator<SAMRecord> query(final String sequence, final int start, final int end, final boolean contained) {
        return mReader.query(sequence, start, end, contained);
    }

    public CloseableIterator<SAMRecord> queryOverlapping(final String sequence, final int start, final int end) {
        return query(sequence, start, end, false);
    }

    public CloseableIterator<SAMRecord> queryContained(final String sequence, final int start, final int end) {
        return query(sequence, start, end, true);
    }

    private void init(final InputStream stream, final boolean eagerDecode) {

        try {
            final BufferedInputStream bufferedStream = toBufferedStream(stream);
            if (isBAMFile(bufferedStream)) {
                mIsBinary = true;
                mReader = new BAMFileReader(bufferedStream, eagerDecode);
            } else if (isSAMFile(bufferedStream)) {
                mIsBinary = false;
                mReader = new SAMTextReader(bufferedStream);
            } else {
                throw new SAMFormatException("Unrecognized file format");
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private void init(final File file, File indexFile, final boolean eagerDecode) {

        try {
            final BufferedInputStream bufferedStream =
                new BufferedInputStream(new FileInputStream(file));
            if (isBAMFile(bufferedStream)) {
                bufferedStream.close();
                mIsBinary = true;
                final BAMFileReader reader = new BAMFileReader(file, eagerDecode);
                mReader = reader;
                if (indexFile == null) {
                    indexFile = findIndexFile(file);
                }
                if (indexFile != null) {
                    mFileIndex = new BAMFileIndex(indexFile);
                    reader.setFileIndex(mFileIndex);
                }
            } else if (isSAMFile(bufferedStream)) {
                if (indexFile != null) {
                    bufferedStream.close();
                    throw new RuntimeException("Cannot use index file with textual SAM file");
                }
                mIsBinary = false;
                mReader = new SAMTextReader(bufferedStream, file);
            } else {
                bufferedStream.close();
                throw new SAMFormatException("Unrecognized file format");
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private File findIndexFile(final File dataFile) {
        final File indexFile =
            new File(dataFile.getParent(), dataFile.getName() + ".bai");
        if (indexFile.exists()) {
            return indexFile;
        } else {
            return null;
        }
    }

    private boolean isBAMFile(final InputStream stream)
        throws IOException {
        return BlockCompressedInputStream.isValidFile(stream);
    }

    private boolean isSAMFile(final InputStream stream) {
        // For now, assume every non-binary file is a SAM text file.
        return true;
    }

    private BufferedInputStream toBufferedStream(final InputStream stream) {
        if (stream instanceof BufferedInputStream) {
            return (BufferedInputStream) stream;
        } else {
            return new BufferedInputStream(stream);
        }
    }
}
