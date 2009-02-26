/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import java.io.File;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;

/**
 * CommandLineProgram to generate to invoke BustardToBamWriter
 *
 * @author Kathleen Tibbetts
 */
public class BustardToSam extends CommandLineProgram {
    // The following attributes define the command-line arguments
    @Usage(programVersion="1.0")
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Generate a BAM binary file from data in an illumina Bustard directory.\n";

    @Option(shortName = "B", doc = "Bustard directory to parse. ")
    public File BUSTARD_DIRECTORY;

    @Option(shortName = "F", doc = "The flowcell. ")
    public String FLOWCELL;

    @Option(shortName = "L", doc = "The lane for which to parse data. ")
    public Integer LANE;

    @Option(shortName = "P", doc = "Whether the lane was a paired-end run. ")
    public Boolean PE;

    @Option(shortName = "O", doc = "The directory for the binary output file. ")
    public File OUTPUT;

    @Override
	protected int doWork() {
        BustardToSamWriter writer = new BustardToSamWriter(
                new BustardFileParser(BUSTARD_DIRECTORY, LANE, PE), OUTPUT, FLOWCELL);
        writer.writeBamFile();
        return 0;
    }

    public static void main(String[] argv) {
        System.exit(new BustardToSam().instanceMain(argv));
    }


}
