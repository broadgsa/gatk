/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.alignment.reference.bwt;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * Writes .ann files - an alternate sequence dictionary format
 * used by BWA/C.  For best results, the input sequence dictionary
 * should be created with Picard's CreateSequenceDictionary.jar,
 * TRUNCATE_NAMES_AT_WHITESPACE=false.
 *
 * @author mhanna
 * @version 0.1
 */
public class ANNWriter {
    /**
     * BWA uses a fixed seed of 11, written into every file.
     */
    private static final int BNS_SEED = 11;

    /**
     * A seemingly unused value that appears in every contig in the ANN.
     */
    private static final int GI = 0;

    /**
     * Input stream from which to read BWT data.
     */
    private final PrintStream out;

    /**
     * Create a new ANNWriter targeting the given file.
     * @param file file into which ANN data should be written.
     * @throws IOException if there is a problem opening the output file.
     */
    public ANNWriter(File file) throws IOException {
        out = new PrintStream(file);
    }

    /**
     * Create a new ANNWriter targeting the given OutputStream.
     * @param stream Stream into which ANN data should be written.
     */
    public ANNWriter(OutputStream stream)  {
        out = new PrintStream(stream);
    }

    /**
     * Write the contents of the given dictionary into the ANN file.
     * Assumes that no ambs (blocks of indeterminate base) are present in the dictionary.
     * @param dictionary Dictionary to write.
     */
    public void write(SAMSequenceDictionary dictionary) {
        long genomeLength = 0L;
        for(SAMSequenceRecord sequence: dictionary.getSequences())
            genomeLength += sequence.getSequenceLength();
        
        int sequences = dictionary.getSequences().size();

        // Write the header
        out.printf("%d %d %d%n",genomeLength,sequences,BNS_SEED);

        for(SAMSequenceRecord sequence: dictionary.getSequences()) {
            String fullSequenceName = sequence.getSequenceName();
            String trimmedSequenceName = fullSequenceName;
            String sequenceComment = "(null)";

            long offset = 0;

            // Separate the sequence name from the sequence comment, based on BWA's definition.
            // BWA's definition appears to accept a zero-length contig name, so mimic that behavior.
            if(fullSequenceName.indexOf(' ') >= 0) {
                trimmedSequenceName = fullSequenceName.substring(0,fullSequenceName.indexOf(' '));
                sequenceComment = fullSequenceName.substring(fullSequenceName.indexOf(' ')+1);
            }

            // Write the sequence GI (?), name, and comment.
            out.printf("%d %s %s%n",GI,trimmedSequenceName,sequenceComment);
            // Write the sequence offset, length, and ambs (currently fixed at 0).
            out.printf("%d %d %d%n",offset,sequence.getSequenceLength(),0);
        }
    }

    /**
     * Close the given output stream.
     */
    public void close() {
        out.close();
    }
}
