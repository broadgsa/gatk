package org.broadinstitute.sting.alignment.reference.bwt;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMSequenceDictionary;

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
