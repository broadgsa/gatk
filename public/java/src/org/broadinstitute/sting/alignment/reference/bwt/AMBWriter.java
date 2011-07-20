package org.broadinstitute.sting.alignment.reference.bwt;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Writes .amb files - a file indicating where 'holes' (indeterminant bases)
 * exist in the contig.  Currently, only empty, placeholder AMBs are supported.
 *
 * @author mhanna
 * @version 0.1
 */
public class AMBWriter {
    /**
     * Number of holes is fixed at zero.
     */
    private static final int NUM_HOLES = 0;

    /**
     * Input stream from which to read BWT data.
     */
    private final PrintStream out;

    /**
     * Create a new ANNWriter targeting the given file.
     * @param file file into which ANN data should be written.
     * @throws java.io.IOException if there is a problem opening the output file.
     */
    public AMBWriter(File file) throws IOException {
        out = new PrintStream(file);
    }

    /**
     * Create a new ANNWriter targeting the given OutputStream.
     * @param stream Stream into which ANN data should be written.
     */
    public AMBWriter(OutputStream stream)  {
        out = new PrintStream(stream);
    }

    /**
     * Write the contents of the given dictionary into the AMB file.
     * Assumes that there are no holes in the dictionary.
     * @param dictionary Dictionary to write.
     */
    public void writeEmpty(SAMSequenceDictionary dictionary) {
        long genomeLength = 0L;
        for(SAMSequenceRecord sequence: dictionary.getSequences())
            genomeLength += sequence.getSequenceLength();

        int sequences = dictionary.getSequences().size();

        // Write the header
        out.printf("%d %d %d%n",genomeLength,sequences,NUM_HOLES);
    }

    /**
     * Close the given output stream.
     */
    public void close() {
        out.close();
    }
}