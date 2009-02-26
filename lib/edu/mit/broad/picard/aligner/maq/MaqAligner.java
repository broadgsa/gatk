/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.aligner.maq;

import edu.mit.broad.picard.aligner.Aligner;
import edu.mit.broad.picard.aligner.AbstractBaseAligner;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.util.Log;

import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

/**
 * Maq implementation of the Aligner interface
 */
public class MaqAligner extends AbstractBaseAligner implements Aligner {

    // Constants related to Maq output files
    public static final String MAQ_MAP_SUFFIX = ".out.aln.map";
    public static final String MAQ_LOG_SUFFIX = ".out.map.log";

    // Internal constant for multi-plexing lane data
    private static final int READ_CHUNK_SIZE = 2000000;

    public static final String REFERENCE_FILE_SUFFIX = ".bfa";

    private final Log log = Log.getInstance(MaqAligner.class);

    private String commandLine = null;
    

    /**
     * Constructor that sets every parameter.  All other constructors delegate to this one.
     *
     * @param stringency            the stringency of the alignment
     * @param readsBamFile          the BAM file containing the reads
     * @param outputPrefix          the directory and filename prefix for output
     * @param referenceFileDir      the directory where the reference file is located
     * @param clipPoints            the clip points
     * @param expectedInsertSize    the expected insert size (null for non-PE lanes)
     * @param readsToAlign          the number of reads to align
     * @param customParametersMap   parameters specific to the Aligner implementation
     */
    public MaqAligner(Stringency stringency, File readsBamFile, String outputPrefix,
                               String referenceFileDir, int clipPoints[], Integer expectedInsertSize,
                               Integer readsToAlign, Map<String, String> customParametersMap,
                               boolean pairedReads, int readLength) {

        super(stringency, readsBamFile, outputPrefix, referenceFileDir, clipPoints,
              expectedInsertSize, readsToAlign, customParametersMap, pairedReads, readLength);
    }

    /**
     * Prepares all the necessary inputs for the alignment process from a BAM file of read data.
     */
    public void prepareInputs() {
        log.info("Preparing Maq inputs.");
        BamToBfqWriter writer = new BamToBfqWriter(this.getReadsBamFile(), this.getOutputPrefix(),
                    this.getReadsToAlign(), READ_CHUNK_SIZE, isPairedReads());
        writer.writeBfqFiles();
    }

    /**
     * Does the alignment and produces output in the underlying form of the aligner.
     */
    public void align() {
        log.info("Running Maq alignment.");

        // Temporary hack until we get the multi-tasking code from Seva
        List<String> mapFileNames = new ArrayList<String>(); // All map files that we will merge together at the end

        String maqParams = MaqConstants.SWITCH_RANDOM_SEED + " " + MaqConstants.DEFAULT_RANDOM_SEED;

        if (this.getStringency() == Stringency.high) {
            maqParams += " " + MaqConstants.SWITCH_MAX_OUTER_DISTANCE + " " + Math.round(
                    this.getExpectedInsertSize() * MaqConstants.HIGH_STRINGENCY_MAX_OUTER_DISTANCE_MULTIPLIER);
            maqParams += " " + MaqConstants.SWITCH_SUM_MISMATCHES + " " +
                    MaqConstants.HIGH_STRINGENCY_SUM_MISMATCHES;
        }
        else {
            maqParams += " " + MaqConstants.SWITCH_MAX_OUTER_DISTANCE + " " +
                    MaqConstants.LOW_STRINGENCY_MAX_OUTER_DISTANCE;
            // For low stringency, get at least 30 bases and then let half of what's remaining mismatch
            int maxMisMatches = (this.getReadLength() - 30)/2;
            maqParams += " " + MaqConstants.SWITCH_SUM_MISMATCHES + " " +
                    (maxMisMatches * MaqConstants.LOW_STRINGENCY_QUALITY_FOR_MISMATCHES);
        }

        String referenceFile = new File(this.getReferenceFileDir()).listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(REFERENCE_FILE_SUFFIX);
            }
        })[0].getAbsolutePath();

        ProcessBuilder builder;

        // Map the bfq files, individually or in pairs
        SortedSet<File> bfqs = new TreeSet<File>(this.getBfqFiles());
        for (Iterator<File> it = bfqs.iterator(); it.hasNext();) {

            String read1bfq = it.next().getAbsolutePath();
            String read2bfq = (this.isPairedReads()) ? it.next().getAbsolutePath() : "";

            String outputFileBase = read1bfq.substring(0, read1bfq.lastIndexOf('.')-2);
            String mapFile = outputFileBase + MAQ_MAP_SUFFIX;
            String logFile = outputFileBase + MAQ_LOG_SUFFIX;

            String command = MaqConstants.MAQ_HOME + MaqConstants.MAQ_COMMAND + " " + MaqConstants.MAP_COMMAND +
                    " " + maqParams + " " + mapFile + " " + referenceFile + " " + read1bfq + " " + read2bfq +
                    " 2> " + logFile;
            setCommandLine(getCommandLine() == null ? command : getCommandLine() + ";" + command);
            log.info("Executing command: " + command);
            try {
                builder = new ProcessBuilder(command.split(" "));
                Process p = builder.start();
                p.waitFor();
            }
            catch (Exception e) {
                throw new PicardException("Error starting Maq process", e);
            }

            mapFileNames.add(mapFile);   
        }

        // If there's more than one map file, then merge them.
        String finalFileName = this.getOutputPrefix() + "." + this.getStringency() + MAQ_MAP_SUFFIX;
        if (mapFileNames.size() > 1) {
            String command = MaqConstants.MAQ_HOME + MaqConstants.MAQ_COMMAND + " " +
                    MaqConstants.MERGE_COMMAND + " " + finalFileName;
            for (String name : mapFileNames) {
                command += " " + name;
            }
            setCommandLine(getCommandLine() == null ? command : getCommandLine() + ";" + command);
            log.info("Executing command: " + command);

            try {
                builder = new ProcessBuilder(command.split(" "));
                Process p = builder.start();
                p.waitFor();
            }
            catch (Exception e) {
                throw new PicardException("Error starting Maq process", e);
            }
        }
        else { // Otherwise rename the single map file so we can find it later
            File f = new File(mapFileNames.get(0));
            if (!f.renameTo(new File(finalFileName))) {
                throw new PicardException("Error renaming " + f.getAbsolutePath() + " to " + finalFileName);
            }
        }
    }

    /**
     * Converts the output of the aligner to BAM format
     */
    public void prepareOutput() {
        log.info("Preparing output from Maq alignment.");
        // TODO: MaqToBam
    }

    /**
     * Cleans up intermediate files (the files created in by and for the underlying aligner by the
     * prepareInputs() and align() methods.  Does not clean up the original source files or the final BAM file.
     */
    public void cleanup() {
        log.info("Cleaning up Maq intermediate files.");
        this.deleteFiles(getBfqFiles());
//        this.deleteFiles(getMaqAlignmentFiles());
    }

    /**
     * Returns a list of zero to two BFQ files, depending on whether they are there
     * and whether it was a paired-end run or not
     *
     * @return a list of BFQ files
     */
    private List<File> getBfqFiles() {
        File dir = new File(this.getOutputPrefix().substring(0, this.getOutputPrefix().lastIndexOf("/")));
        return Arrays.asList(dir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(".bfq");
            }
        }));
    }

    /**
     * Returns the Maq map files
     *
     * @return a list of Maq .map files
     */
    private List<File> getMaqAlignmentFiles() {
        File dir = new File(this.getOutputPrefix().substring(0, this.getOutputPrefix().lastIndexOf("/")));
        return Arrays.asList(dir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                // TODO: Add the text files if we do not read the binary map files
                return name.endsWith(MAQ_MAP_SUFFIX) || name.endsWith(MAQ_LOG_SUFFIX);
            }
        }));
    }

    public String getCommandLine() { return commandLine; }
    public void setCommandLine(String commandLine) { this.commandLine = commandLine; }
}
