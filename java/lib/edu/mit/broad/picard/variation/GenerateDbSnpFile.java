/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.variation;

import java.io.File;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;

/**
 * CommandLineProgram to generate to invoke DbSnpFileGenerator
 *
 * @author Kathleen Tibbetts 
 */
public class GenerateDbSnpFile extends CommandLineProgram
{
    // The following attributes define the command-line arguments
    @Usage(programVersion="1.0")
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Generate a KnownVariant binary file from a UCSC DbSnp text file.\n";

    @Option(shortName = "S", doc = "UCSC SNP file. ")
    public File SNP_FILE;

    @Option(shortName = "D", doc = "Sequence Dictionary for the genome in SAM or BAM format. ")
    public File SEQUENCE_DICTIONARY;

    @Option(shortName = "O", doc = "The binary output file. ")
    public File OUTPUT;

    @Override
	protected int doWork() {
        DbSnpFileGenerator generator = new DbSnpFileGenerator(SNP_FILE, SEQUENCE_DICTIONARY, OUTPUT);
        generator.writeDbSnpFile();
        return 0; 
    }

    public static void main(String[] argv) {
        System.exit(new GenerateDbSnpFile().instanceMain(argv));
    }

}
