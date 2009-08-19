package org.broadinstitute.sting.utils.genotype.vcf;


import java.io.File;


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
        if (args.length != 1) {
            printUsage();
            return;
        }
        File vcfFile = new File(args[0]);
        if (!vcfFile.exists()) {
            System.err.println("Specified VCF file doesn't exist, please check the input file\n");
            printUsage();
            return;
        }
        int counter = 0;
        try {
            VCFReader reader = new VCFReader(vcfFile);
            while (reader.hasNext()) {
                counter++;
                reader.next();
            }
        } catch (Exception e) {
            System.err.println("VCF Validation failed, after parsing " + counter + " entries.");
            System.err.println("The reason given was: " + e.getMessage());            
        }
        System.err.println("Viewed " + counter + " VCF record entries.");
    }

    public static void printUsage() {
        System.err.println("VCF validator (VCF Version " + VCF_VERSION + ")");
        System.err.println("Usage:");
        System.err.println("vcfvalidator <fille.vcf>");
        System.err.println("");
    }

}
