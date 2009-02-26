package edu.mit.broad.picard.reference;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMTextHeaderCodec;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.util.LineReader;
import edu.mit.broad.sam.util.AsciiLineReader;

import java.io.*;
import java.nio.charset.Charset;
import java.util.List;

/**
 * Implementation of ReferenceSequenceFile for reading from FASTA files.
 *
 * @author Tim Fennell
 */
class FastaSequenceFile implements ReferenceSequenceFile {
    private static final Charset ASCII = Charset.forName("US-ASCII");
    private File file;
    private BufferedReader in;
    private List<SAMSequenceRecord> sequenceDictionary;
    private String cachedLine = null;
    private int index = -1;

    /** Constructs a FastaSequenceFile that reads from the specified file. */
    FastaSequenceFile(File file) {
        this.file = file;
        this.in = new BufferedReader(new InputStreamReader(IoUtil.openFileForReading(file)));

        // Try and locate the dictionary
        String dictionaryName = file.getAbsolutePath();
        dictionaryName = dictionaryName.substring(0, dictionaryName.lastIndexOf(".fasta"));
        dictionaryName += ".dict";
        File dictionary = new File(dictionaryName);
        if (dictionary.exists()) {
            IoUtil.assertFileIsReadable(dictionary);

            try {
                SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                SAMFileHeader header = codec.decode(new AsciiLineReader(new FileInputStream(dictionary)), dictionary);
                if (header.getSequences() != null && header.getSequences().size() > 0) {
                    this.sequenceDictionary = header.getSequences();
                }
            }
            catch (Exception e) {
                throw new PicardException("Could not open sequence dictionary file: " + dictionaryName, e);
            }
        }
    }

    /**
     * Returns the list of sequence records associated with the reference sequence if found
     * otherwise null.
     */
    public List<SAMSequenceRecord> getSequenceDictionary() {
        return this.sequenceDictionary;
    }

    public ReferenceSequence nextSequence() {
        String line = null;
        String name = null;

        // Scan forward to a header line
        while ((line = readNextLine()) != null) {
            if (line.startsWith(">")) {
                name = line.substring(1).trim();
                this.index += 1;
                break;
            }
        }

        // No more!
        if (name == null) return null;

        // Read the sequence
        int basesRead = 0;
        byte[] bases = new byte[250000000]; // big enough to hold human chr1!
        while ((line = readNextLine()) != null) {
            if (line.startsWith(">")) {
                pushBackLine(line);
                break;
            }
            else {
                final byte[] nextBases = line.getBytes(ASCII);
                final int lineLength = nextBases.length;

                // If the array isn't big enough to hold the next chunk, resize it
                if (basesRead + lineLength > bases.length) {
                    byte[] tmp = new byte[bases.length * 2];
                    System.arraycopy(bases, 0, tmp, 0, basesRead);
                    bases = tmp;
                }

                // Now shunt the most recent bases onto the end of the array
                System.arraycopy(nextBases, 0, bases, basesRead, lineLength);
                basesRead += lineLength;
            }
        }

        // And lastly resize the array down to the right size
        if (basesRead != bases.length) {
            byte[] tmp = new byte[basesRead];
            System.arraycopy(bases, 0, tmp, 0, basesRead);
            bases = tmp;
        }

        return new ReferenceSequence(name, this.index, bases);
    }

    /**
     * Reads the next line from the file, or if we've saved a line earlier, returns that
     * instead.
     */
    private String readNextLine() {
        // If we have a cached line use it
        if (this.cachedLine != null) {
            String tmp = this.cachedLine;
            this.cachedLine = null;
            return tmp;
        }
        else {
            try { return this.in.readLine(); }
            catch (IOException ioe) {
                throw new PicardException("Error reading line from file: " + this.file.getAbsolutePath(), ioe);
            }
        }
    }

    /** Pushed a line back so that the next call to readNextLine() will return it. */
    private void pushBackLine(String line) {
        this.cachedLine = line;
    }
}

