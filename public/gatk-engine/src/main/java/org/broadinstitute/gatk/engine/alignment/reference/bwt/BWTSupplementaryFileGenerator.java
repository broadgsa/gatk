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

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.SAMSequenceDictionary;

import java.io.File;
import java.io.IOException;

/**
 * Generate BWA supplementary files (.ann, .amb) from the command line.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWTSupplementaryFileGenerator {
    enum SupplementaryFileType { ANN, AMB } 

    public static void main(String[] args) throws IOException {
        if(args.length < 3)
            usage("Incorrect number of arguments supplied");

        File fastaFile = new File(args[0]);
        File outputFile = new File(args[1]);
        SupplementaryFileType outputType = null;
        try {
            outputType = Enum.valueOf(SupplementaryFileType.class,args[2]);
        }
        catch(IllegalArgumentException ex) {
            usage("Invalid output type: " + args[2]);
        }

        ReferenceSequenceFile sequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile);
        SAMSequenceDictionary dictionary = sequenceFile.getSequenceDictionary();

        switch(outputType) {
            case ANN:
                ANNWriter annWriter = new ANNWriter(outputFile);
                annWriter.write(dictionary);
                annWriter.close();
                break;
            case AMB:
                AMBWriter ambWriter = new AMBWriter(outputFile);
                ambWriter.writeEmpty(dictionary);
                ambWriter.close();
                break;
            default:
                usage("Unsupported output type: " + outputType);
        }
    }

    /**
     * Print usage information and exit.
     */
    private static void usage(String message) {
        System.err.println(message);
        System.err.println("Usage: BWTSupplementaryFileGenerator <fasta> <output file> <output type>");
        System.exit(1);
    }
}
