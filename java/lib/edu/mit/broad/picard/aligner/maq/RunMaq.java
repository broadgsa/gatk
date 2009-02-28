/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner.maq;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.aligner.Aligner;

import java.io.File;
import java.util.Map;
import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;

/**
 * CommandLineProgram to generate to invoke BustardToBamWriter
 *
 * @author Kathleen Tibbetts
 */
public class RunMaq extends CommandLineProgram {
    private static final String PROGRAM_VERSION = "1.0";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Invoke the Maq aligner.\n" +
            "Version: " + PROGRAM_VERSION +"\n";

    @Option(shortName="I", doc="The BAM file to parse.", optional=true)
    public File INPUT;
    @Option(shortName="O", doc="The directory and file prefix for all output.", optional=false)
    public String OUTPUT;
    @Option(shortName="L", doc="The read length.", optional=false)
    public Integer READ_LENGTH;
    @Option(shortName="S", doc="Stringency of the alignment.", optional=true)
    public Aligner.Stringency STRINGENCY;
    @Option(shortName="R", doc="Directory where the reference file is located.", optional=true)
    public String REFERENCE;
    @Option(shortName="C", doc="Clip points for the alignment.", optional=true, minElements=0, maxElements=4)
    public List<Integer> CLIP_POINT = new ArrayList<Integer>();
    @Option(shortName="E", doc="Expected insert size.", optional=true)
    public Integer EXPECTED_INSERT_SIZE;
    @Option(doc="Whether this is a paired-end run.", optional=false)
    public Boolean PE;
    @Option(shortName="NUM", doc="Number of reads to align (null = all).", optional=true)
    public Integer READS_TO_ALIGN;
    @Option(shortName="CUSTOM", doc="Custom parameter in the form name=value.", optional=true)
    public List<String> CUSTOM_PARAMETER = new ArrayList<String>();
    @Option(shortName="PREP", doc="Whether to prepare inputs for the alignement.", optional=true)
    public Boolean PREPARE = true;
    @Option(doc="Whether to do the alignement.", optional=true)
    public Boolean ALIGN = true;
    @Option(shortName="BAM", doc="Whether to generate a BAM file from the alignment output.", optional=true)
    public Boolean BAM_OUTPUT = true;
    @Option(doc="Whether to clean up intermediate input and output.", optional=true)
    public Boolean CLEANUP = true;

    protected int doWork() {
        int clipPoints[] = null;
        if (CLIP_POINT != null) {
            clipPoints = new int[4];
            int index=0;
            for (Integer i : CLIP_POINT) {
                clipPoints[index++] = i;
            }
        }
        Map<String, String> params = null;
        if (CUSTOM_PARAMETER != null) {
            params = new HashMap<String, String>();
            for (String param : CUSTOM_PARAMETER) {
                String nameAndVal[] = param.split("=");
                params.put(nameAndVal[0], nameAndVal[1]);
            }
        }
        Aligner aligner = new MaqAligner(STRINGENCY, INPUT, OUTPUT, REFERENCE, clipPoints,
                EXPECTED_INSERT_SIZE, READS_TO_ALIGN, params, PE, READ_LENGTH);
        if (PREPARE) {
            aligner.prepareInputs();
        }
        if (ALIGN) {
            aligner.align();
        }
        if (BAM_OUTPUT) {
            aligner.prepareOutput();
        }
        if (CLEANUP) {
            aligner.cleanup();
        }
        return 0;
    }

    /**
     * This is kind of a mess.  Almost everything is optional, since you don't have to do all of the steps in the
     * alignement.
     * @return
     */
    protected boolean customCommandLineValidation() {
        if (PREPARE) {
            if( INPUT == null) {
                System.err.println("ERROR: INPUT must be specified when preparing inputs for the alignment.");
                return false;
            }
            if (CLIP_POINT.size() != 0 && CLIP_POINT.size() != 4) {
                System.err.println("ERROR: You must supply either 0 or 4 values for CLIP_POINT: " + CLIP_POINT.size());
                return false;
            }
        }
        if (ALIGN) {
            if (STRINGENCY == null) {
                System.err.println("ERROR: STRINGENCY must be specified when doing an alignment.");
                return false;
            }
            if (REFERENCE == null) {
                System.err.println("ERROR: REFERENCE must be specified when doing an alignment.");
                return false;
            }
        }
        return true;
    }

    public static void main(String[] argv) {
        System.exit(new RunMaq().instanceMain(argv));
    }
}
