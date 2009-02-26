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


import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.util.BlockCompressedInputStream;
import edu.mit.broad.sam.util.CloseableIterator;
import edu.mit.broad.sam.util.RuntimeIOException;


/**
 * Class for reading GELI (GEnotype LIkelihood) files.
 * 
 * @author Doug Voet
 */
public class GeliFileReader implements Iterable<GenotypeLikelihoods>
{
    private ReaderImplementation mReader = null;

    /**
     * Internal interface for SAM/BAM file reader implementations.
     * Implemented as an abstract class to enforce better access control.
     */
    static abstract class ReaderImplementation {
        abstract SAMFileHeader getFileHeader();
        abstract CloseableIterator<GenotypeLikelihoods> getIterator();
        abstract void close();
    }


    public GeliFileReader(final InputStream stream) {
        try {
            final BufferedInputStream bufferedStream = toBufferedStream(stream);
            if (isValidGELIFile(bufferedStream)) {
                mReader = new GeliFileReaderImplementation(bufferedStream);
            } else {
                throw new GeliException("Unrecognized file format");
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    public GeliFileReader(final File file) {
        try {
            final BufferedInputStream bufferedStream =
                new BufferedInputStream(new FileInputStream(file));
            if (isValidGELIFile(bufferedStream)) {
                bufferedStream.close();
                final GeliFileReaderImplementation reader = new GeliFileReaderImplementation(file);
                mReader = reader;
            } else {
                bufferedStream.close();
                throw new GeliException("Unrecognized file format");
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    public void close() {
        if (mReader != null) {
            mReader.close();
        }
        mReader = null;
    }

    public SAMFileHeader getFileHeader() {
        return mReader.getFileHeader();
    }

    public CloseableIterator<GenotypeLikelihoods> iterator() {
        return mReader.getIterator();
    }

    private boolean isValidGELIFile(final InputStream stream)
        throws IOException {
        return BlockCompressedInputStream.isValidFile(stream);
    }

    private BufferedInputStream toBufferedStream(final InputStream stream) {
        if (stream instanceof BufferedInputStream) {
            return (BufferedInputStream) stream;
        } else {
            return new BufferedInputStream(stream);
        }
    }
}
