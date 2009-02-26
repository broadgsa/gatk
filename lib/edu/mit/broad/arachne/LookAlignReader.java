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
package edu.mit.broad.arachne;


import edu.mit.broad.sam.util.CloseableIterator;

import java.io.*;


/**
 * Reader for arachne LookAlign text format alignment files.
 * Supports filtering of the input by genomic locus.
 */
public class LookAlignReader
    implements CloseableIterator<Alignment> {

    private LineNumberReader mReader = null;
    private Alignment mNextAlignment = null;
    private int mBSequenceId = -1;
    private int mBStart = 0;
    private int mBEnd = 0;


    public LookAlignReader(File file)
        throws IOException {
        this(new FileReader(file));
    }

    public LookAlignReader(Reader reader) {
        if (reader instanceof LineNumberReader) {
            mReader = (LineNumberReader) reader;
        } else {
            mReader = new LineNumberReader(reader);
        }
    }

    public void setBSequenceId(int value) {
        mBSequenceId = value;
    }

    public void setBStart(int value) {
        mBStart = value;
    }

    public void setBEnd(int value) {
        mBEnd = value;
    }

    public boolean hasNext() {
        if (mNextAlignment != null) {
            return true;
        }
        try {
            mNextAlignment = nextAlignment();
            return (mNextAlignment != null);
        } catch (IOException exc) {
            throw new RuntimeException(exc.getMessage(), exc);
        }
    }

    public Alignment next() {
        if (!hasNext()) {
            throw new IllegalStateException("Iterator exhausted");
        }
        try {
            Alignment result = mNextAlignment;
            mNextAlignment = nextAlignment();
            return result;
        } catch (IOException exc) {
            throw new RuntimeException(exc.getMessage(), exc);
        }
    }

    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }

    public void close() {
        if (mReader != null) {
            try {
                mReader.close();
            } catch (IOException exc) {
                throw new RuntimeException(exc.getMessage(), exc);
            }
            mReader = null;
        }
    }

    private Alignment nextAlignment()
        throws IOException {
        if (mReader == null) {
            return null;
        }
        while (true) {
            String line = mReader.readLine();
            if (line == null) {
                close();
                break;
            }
            if (!line.startsWith("QUERY")) {
                continue;
            }
            Alignment alignment = Alignment.parse(line);
            if (matchesFilters(alignment)) {
                return alignment;
            }
        }
        return null;
    }

    private boolean matchesFilters(Alignment alignment) {
        if (mBSequenceId < 0) {
            return true;
        }
        if (alignment.getBSequenceId() != mBSequenceId) {
            return false;
        }
        if (mBStart > 0 && alignment.getBEnd() < mBStart) {
            return false;
        }
        if (mBEnd > 0 && alignment.getBStart() > mBEnd) {
            return false;
        }
        return true;
    }
}

