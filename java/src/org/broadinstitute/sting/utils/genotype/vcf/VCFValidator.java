package org.broadinstitute.sting.utils.genotype.vcf;


import java.io.File;
import java.util.Map;
import java.util.TreeMap;


/**
 * @author aaron
 *         <p/>
 *         Class VCFValidator
 *         <p/>
 *         This is the main class for providing a light weight validation of a VCF file.
 *         It has two parameters, an optional -A flag meaning that you'd like to collect all
 *         the errors and present them at the end, and the VCF file itself (a required parameter).
 */
public class VCFValidator {

    private static final String VCF_VERSION = "VCFv3.2";

    /**
     * about as simple as things come right now.  We open the file, process all the entries in the file,
     * and if no errors pop up in processing, well hey, looks good to us.
     *
     * @param args the vcf file is the only required parameter, with the optional -A indicating that errors
     * should be held until the end of processing
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
        Map<Integer, Exception> problems = new TreeMap<Integer, Exception>();

        try {
            // open up our reader
            VCFReader reader = new VCFReader(vcfFile);

            // the number of samples should be set in the header and consistant over all records
            final int sampleCount = reader.getHeader().getGenotypeSamples().size();
            while (reader.hasNext()) {
                try {
                    recordCount++;
                    VCFRecord rec = reader.next();
                    // if the header indicates we have genotyping data, try to extract it for all samples
                    if (reader.getHeader().hasGenotypingData()) {
                        int sampleCounter = 0;
                        for (VCFGenotypeRecord genorec : rec.getVCFGenotypeRecords()) {
                            sampleCounter++;
                            /**
                             * just cycle through the records right now; any additional checks for
                             * the records should go in this block.
                             **/
                        }
                        if (sampleCounter != sampleCount)
                            throw new RuntimeException("Record " + recordCount + " does not have the required number " +
                                    "of records (" + sampleCounter + " in the record, " + sampleCount + " in the header)");
                        
                    }
                } catch (Exception e) {
                    if (catchAll)
                        problems.put(recordCount, e);
                    else {
                        validationFailed(e, recordCount);
                        return;
                    }
                }
            }
        } catch (Exception e) {
            if (catchAll)
                problems.put(new Integer(0), e);
            else
                validationFailed(e, recordCount);
        }
        System.err.println("Viewed " + recordCount + " VCF record entries.");
        if (problems.size() > 0) {
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

    /**
     * print the usage information for the VCF validator
     */
    public static void printUsage() {
        System.err.println("VCF validator (VCF Version " + VCF_VERSION + ")");
        System.err.println("Usage:");
        System.err.println("vcfvalidator <-A> <fille.vcf>");
        System.err.println("");
        System.err.println("\t-A\tTell the validator to attempt to catch all the problems, and not stop at the first.  Some may be too fatal to continue.");
        System.err.println("");
    }

}
