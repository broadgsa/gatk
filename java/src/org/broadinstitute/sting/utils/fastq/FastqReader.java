package org.broadinstitute.sting.utils.fastq;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.exceptions.UserError;

import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.io.*;

public class FastqReader implements Iterator<FastqRecord>, Iterable<FastqRecord>, Closeable {
    private File fastqFile;
    private BufferedReader in;
    private FastqRecord nextRecord;

    public FastqReader(File file) {
        fastqFile = file;

        try {
            if (fastqFile.getName().endsWith(".gz")) {
                in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastqFile))));
            } else {
                in = new BufferedReader(new InputStreamReader(new FileInputStream(fastqFile)));
            }

            nextRecord = readNextRecord();
        } catch (IOException e) {
            throw new UserError.CouldNotReadInputFile(file, String.format("Error opening '%s'", fastqFile.getAbsolutePath()));
        }
    }

    private FastqRecord readNextRecord() {
        try {
            String seqHeader = in.readLine();
            String seqLine = in.readLine();
            String qualHeader = in.readLine();
            String qualLine = in.readLine();

            return new FastqRecord(seqHeader, seqLine, qualHeader, qualLine);
        } catch (IOException e) {
            throw new UserError.CouldNotReadInputFile(fastqFile, e);
        }
    }

    public boolean hasNext() { return nextRecord != null; }

    public FastqRecord next() {
        FastqRecord rec = nextRecord;

        try {
            if (in.ready()) {
                nextRecord = readNextRecord();
            } else {
                nextRecord = null;
            }
        } catch (IOException e) {
            throw new UserError.CouldNotReadInputFile(fastqFile, e);
        }

        return rec;
    }

    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }

    public Iterator<FastqRecord> iterator() { return this; }

    public void close() {
        try {
            in.close();
        } catch (IOException e) {
            throw new UserError.CouldNotReadInputFile(fastqFile, e);
        }
    }
}
