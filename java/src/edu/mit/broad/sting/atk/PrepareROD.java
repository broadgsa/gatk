package edu.mit.broad.sting.atk;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMSequenceRecord;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;

import edu.mit.broad.sting.atk.modules.*;
import edu.mit.broad.sting.utils.*;

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

    public static HashMap<String, Class> Types = new HashMap<String,Class>();
    public static void addModule(final String name, final Class rodType) {
        System.out.printf("* Adding rod class %s%n", name);
        Types.put(name.toLowerCase(), rodType);
    }

    static {
        addModule("GFF", rodGFF.class);
        addModule("dbSNP", rodDbSNP.class);
    }

    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new PrepareROD().instanceMain(argv));
    }

    protected int doWork() {

        // Prepare the sort ordering w.r.t. the sequence dictionary
        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE_ARG);
        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary();
        HashMap<String, Integer> refContigOrdering = new HashMap<String, Integer>();

        int i = 0;
        for ( SAMSequenceRecord contig : refContigs ) {
            System.out.println(contig.getSequenceName());
            refContigOrdering.put(contig.getSequenceName(), i);
            i++;
        }
        GenomeLoc.setContigOrdering(refContigOrdering);

        Class rodClass = Types.get(ROD_TYPE.toLowerCase());

        ReferenceOrderedData rod = new ReferenceOrderedData(new File(ROD_FILE), rodClass );
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
        ReferenceOrderedData outputRod = new ReferenceOrderedData(new File(OUTPUT_FILE), rodClass );
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
