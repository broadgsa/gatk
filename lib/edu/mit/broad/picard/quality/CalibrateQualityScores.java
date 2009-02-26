package edu.mit.broad.picard.quality;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.variation.DbSnpFileReader;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMFileWriter;
import edu.mit.broad.sam.SAMFileWriterFactory;
import edu.mit.broad.sam.SAMRecord;

import java.io.File;
import java.io.PrintStream;

/**
 * Command line program to calibrate quality scores using alignment and dbsnp data. Calibrates
 * qualities cycle by cycle and separately for reads one and two in a pair. Bases that fall
 * within dbSNP loci are ignored otherwise the empircal mismatch rate is calculated for
 * each quality at each cycle and used to calculate the calibrated quality value.
 *
 * @author Tim Fennell
 */
public class CalibrateQualityScores extends CommandLineProgram {
    @Option(shortName="A", doc="A file of aligned reads in SAM or BAM format")
    public File ALIGNED_SAM;

    @Option(shortName="I", doc="A SAM or BAM file to rewrite with calibrated qualities. If omitted ALIGNED_SAM is used.", optional=true)
    public File INPUT;

    @Option(shortName="O", doc="The SAM or BAM file to write with updated qualities.")
    public File OUTPUT;

    @Option(shortName="R", doc="Reference sequence file")
    public File REFERENCE;

    @Option(shortName="SNP", doc="Binary file of dbSNP information", optional=true)
    public File DBSNP_FILE;

    @Option(shortName="TABLE", doc="A file to output the calibration table(s) to.")
    public File CALIBRATION_TABLE_OUT;

    @Option(doc="Optional limit to the number of aligned reads that should be procesed", optional=true)
    public Integer READ_LIMIT = -1;

    /** Stock main method for a command line program. */
    public static void main(String[] argv) {
        System.exit(new CalibrateQualityScores().instanceMain(argv));
    }

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then loads up all the
     * data and calibrates the quality scores and proceeds to write an output file
     * with calibrated quality scores instead of the input quality scores.
     */
    protected int doWork() {
        final Log log = Log.getInstance(getClass());

        // Some quick parameter checking
        if (INPUT == null) INPUT = ALIGNED_SAM;

        IoUtil.assertFileIsReadable(ALIGNED_SAM);
        IoUtil.assertFileIsReadable(REFERENCE);
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(CALIBRATION_TABLE_OUT);

        log.info("Reading input files and calculating calibration matrices.");

        // Load things up and calculate the quality score calibrations
        SAMFileReader sam = new SAMFileReader(ALIGNED_SAM);
        ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        DbSnpFileReader dbsnp = null;

        if (DBSNP_FILE != null) {
            IoUtil.assertFileIsReadable(DBSNP_FILE);
            dbsnp = new DbSnpFileReader(DBSNP_FILE);
        }

        QualityScoreCalibrator calibrator = new QualityScoreCalibrator(sam, ref, dbsnp);
        calibrator.calibrate(READ_LIMIT);

        // Dump the calibration tables
        log.info("Writing out calibration table.");
        PrintStream stream = new PrintStream(IoUtil.openFileForWriting(CALIBRATION_TABLE_OUT));
        stream.println("Read 1 Calibration Table:");
        print(stream, calibrator.getRead1Matrix().getCalibratedQualities());

        if (!calibrator.getRead2Matrix().isEmpty()) {
            stream.println();
            stream.println("Read 2 Calibration Table:");
            print(stream, calibrator.getRead2Matrix().getCalibratedQualities());
        }

        // And then load up the input and rewrite with calibrated qualities
        log.info("Writing file with calibrated qualities.");
        SAMFileReader in  = new SAMFileReader(INPUT);
        SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);

        for (SAMRecord rec : in) {
            byte[] quals = rec.getBaseQualities();
            byte[] calibrated = new byte[quals.length];
            QualityScoreMatrix matrix = rec.getFirstOfPairFlag() ? calibrator.getRead1Matrix() : calibrator.getRead2Matrix();

            for (int i=0; i<quals.length; ++i) {
                calibrated[i] = (byte) matrix.getCalibratedQuality(i+1, quals[i]);
            }

            rec.setBaseQualities(calibrated);
            out.addAlignment(rec);
        }

        out.close();

        return 0;
    }

    /** Static helper method to dump a calibration matrix to the screen for debugging. */
    private void print(PrintStream out, int[][] matrix) {
        int maxY = 0;
        for (int x=0; x<matrix.length; ++x) {
            if (matrix[x] != null) {
                maxY = Math.max(maxY, matrix[x].length - 1);
            }
        }

        // Print out the header row
        for (int i=0;i<=maxY; ++i) {
            out.print(i + "\t");
        }
        out.println();

        // Now print out the data cycle by cycle
        for (int cycle=1; cycle<matrix.length; ++cycle) {
            out.print(cycle + "\t");

            int[] quals = matrix[cycle];

            for (int qual=1; qual<quals.length; ++qual) {
                out.print(quals[qual] + "\t");
            }
            out.println();
        }
    }
}
