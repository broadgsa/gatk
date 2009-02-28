package edu.mit.broad.picard.util;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.io.IoUtil;

import java.io.File;
import java.util.Date;
import java.text.SimpleDateFormat;

/**
 * CommandLineProgram to create Picard analysis directory
 *
 * @author Kathleen Tibbetts
 */
public class CreateAnalysisDirectory extends CommandLineProgram {

    public static final SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd");

    // The following attributes define the command-line arguments
    @Usage(programVersion="1.0")
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Create a new Picard analysis directory.\n";


    @Option(shortName = "P", doc = "Analysis directory prefix.  ")
    public String PREFIX = "/seq/picard";

    @Option(shortName = "F", doc = "The flowcell.  ")
    public String FLOWCELL;

    @Option(shortName = "A", doc = "The first cycle being analyzed.  ")
    public Integer FIRST_CYCLE = 1;

    @Option(shortName = "O", doc = "The last cycle being analyzed.  ")
    public Integer LAST_CYCLE;

    @Option(shortName = "R", doc = "The run date in the format MM/dd/yyyy. ")
    public Date RUNDATE;

    @Option(shortName = "L", doc = "Lane number. ")
    public Integer LANE;

    @Option(shortName="LIB", doc = "Library this analysis is for (e.g. 'Solexa-1234').  ")
    public String LIBRARY;

    @Option(shortName="S", doc = "Analysis start date in the format MM/dd/yyyy")
    public Date ANALYSIS_START_DATE;

    @Override
	protected int doWork() {
        if (PREFIX.charAt(PREFIX.length()-1) == '/') {
            PREFIX = PREFIX.substring(0, PREFIX.length()-1);
        }
        IoUtil.assertDirectoryIsWritable(new File(PREFIX));
        String parts[] = { PREFIX, FLOWCELL, "C" + FIRST_CYCLE + "-" + LAST_CYCLE + "_" +
                dateFormat.format(RUNDATE) + "_" + dateFormat.format(ANALYSIS_START_DATE),
                String.valueOf(LANE), LIBRARY };
        String directory = null;

        for (int i = 1; i < parts.length; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j=0; j <= i; j++) {
                sb.append(parts[j]).append("/");
            }
            directory = sb.toString();
            File dir = new File(directory);
            if (!dir.exists()) {
                if (!dir.mkdir()) {
                    System.err.println("Unable to create directory " + directory);
                    return 1;
                }
            }
        }
        System.out.print(directory);
        return 0;
    }

    public static void main(String[] argv) {
        CreateAnalysisDirectory cmd = new CreateAnalysisDirectory();
        cmd.QUIET = true;
        System.exit(cmd.instanceMain(argv));
    }

    
}
