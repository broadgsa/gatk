package org.broadinstitute.sting.utils.genotype.vcf;


import java.io.File;
import java.util.Map;
import java.util.TreeMap;


/**
 * @author aaron
 *         <p/>
 *         Class VCFValidator
 *         <p/>
 *         validate a VCF file
 */
public class VCFValidator {

    private static final String VCF_VERSION = "VCFv3.2";

    /**
     * about as simple as things come right now.  We open the file, process all the entries in the file,
     * and if no errors pop up in processing, well hey, looks good to us.
     * TODO: add validation to individual records fields as they make sense
     *
     * @param args the vcf file is the only parameter
     */
    public static void main(String[] args) {
        boolean catchAll = false;

        if (args.length == 2 && args[0].equals("-A"))
            catchAll = true;
        else if (args.length == 1)
            catchAll = false;
        else {
            printUsage();
            return;
        }
        File vcfFile = new File(args[(catchAll) ? 1 : 0]);
        if (!vcfFile.exists()) {
            System.err.println("Specified VCF file doesn't exist, please check the input file\n");
            printUsage();
            return;
        }
        // count hom many records we see
        int recordCount = 0;
        Map<Integer,Exception> problems = new TreeMap<Integer,Exception>();

        try {
            // open up our reader
            VCFReader reader = new VCFReader(vcfFile);

            while (reader.hasNext()) {
                try {
                    recordCount++;
                    VCFRecord rec = reader.next();
                    // if the header indicates we have genotyping data, try to extract it for all samples
                    if (reader.getHeader().hasGenotypingData()) {
                        for (VCFGenotypeRecord genorec : rec.getVCFGenotypeRecords()) {
                            // just cycle through them, more checks go here
                        }
                    }
                } catch (Exception e) {
                    if (catchAll)
                        problems.put(recordCount,e);
                    else {
                        validationFailed(e, recordCount);
                        return;
                    }
                }
            }
        } catch (Exception e) {
            if (catchAll)
                problems.put(new Integer(0),e);
            else
                validationFailed(e, recordCount);
        }
        System.err.println("Viewed " + recordCount + " VCF record entries.");
        if (problems.size() > 0)  {
            System.err.println("Encountered " + problems.size() + " number of issues. (record zero indicates a header problem)");
            for (Integer e : problems.keySet()) {
                System.err.println("\tProblem at record " + e + " : " + problems.get(e));               
            }
        }
    }

    /**
     * validation failed
     *
     * @param e     the exception
     * @param count the current record count
     */
    public static void validationFailed(Exception e, int count) {
        System.err.println("VCF Validation failed, after parsing " + count + " entries.");
        System.err.println("The reason given was: " + e.getMessage());
        e.printStackTrace();
    }

    /** print the usage information for the VCF validator */
    public static void printUsage() {
        System.err.println("VCF validator (VCF Version " + VCF_VERSION + ")");
        System.err.println("Usage:");
        System.err.println("vcfvalidator <-A> <fille.vcf>");
        System.err.println("");
        System.err.println("\t-A\tTell the validator to attempt to catch all the problems, and not stop at the first.  Some may be too fatal to continue.");
        System.err.println("");
    }

}
