package org.broadinstitute.sting.playground.gatk;

import net.sf.samtools.SAMSequenceRecord;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequenceFile;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.refdata.*;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

public class PrepareROD extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="0.1") public String USAGE = "SAM Validator\n";
    @Option(shortName="REF", doc="Reference sequence file") public File REF_FILE_ARG = null;
    @Option(shortName="ROD", doc="Referenced Ordered Data file") public String ROD_FILE = null;
    @Option(shortName="OUT", doc="Referenced Ordered Data file") public String OUTPUT_FILE = null;
    @Option(shortName="RODNAME", doc="Name of the data") public String ROD_NAME = null;
    @Option(shortName="RODTYPE", doc="Referenced Ordered Data type") public String ROD_TYPE = null;

    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new PrepareROD().instanceMain(argv));
    }

    protected int doWork() {

        // Prepare the sort ordering w.r.t. the sequence dictionary
        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE_ARG);
        GenomeLoc.setupRefContigOrdering(refFile);

        Class<? extends ReferenceOrderedDatum> rodClass = ReferenceOrderedData.Types.get(ROD_TYPE.toLowerCase()).type;

        ReferenceOrderedData<? extends ReferenceOrderedDatum> rod = new ReferenceOrderedData("ROD", new File(ROD_FILE), rodClass );
        try {
            rod.validateFile();
        } catch ( Exception e ) {
            //System.out.println("Validation failure: " + e);
            e.printStackTrace();
        }

        ArrayList<ReferenceOrderedDatum> rodData = rod.readAll();
        System.out.printf("Read %d elements from %s%n", rodData.size(), ROD_FILE);
        ReferenceOrderedData.sortRODDataInMemory(rodData);
        try {
            ReferenceOrderedData.write(rodData, new File(OUTPUT_FILE));
         } catch ( IOException e ) {
            //System.out.println("Validation failure: " + e);
            e.printStackTrace();
        }

        System.out.printf("Validating output file %s%n", rodData.size(), OUTPUT_FILE);
        ReferenceOrderedData outputRod = new ReferenceOrderedData("ROD", new File(OUTPUT_FILE), rodClass );
        try {
            outputRod.validateFile();
            //outputRod.hasSameContents(ROD_FILE);
        } catch ( Exception e ) {
            //System.out.println("Validation failure: " + e);
            e.printStackTrace();
        }        

        return 0;
    }
}
