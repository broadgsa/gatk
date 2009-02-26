package edu.mit.broad.picard.directed;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.metrics.MetricsFile;
import edu.mit.broad.sam.SAMFileReader;

import java.io.File;

/**
 * Calculates a set of HS metrics from a sam or bam file.
 *
 * @author Tim Fennell
 */
public class CalculateHsMetrics extends CommandLineProgram {
    @Usage public final String USAGE =
            "Calculates a set of Hybrid Selection specific metrics from an aligned SAM" +
            "or BAM file.";
    @Option(shortName="BI") public File BAIT_INTERVALS;
    @Option(shortName="TI") public File TARGET_INTERVALS;
    @Option(shortName="I") public File INPUT;
    @Option(shortName="M") public File METRICS_FILE;

    /** Stock main method. */
    public static void main(String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    /**
     * Asserts that files are readable and writable and then fires off an
     * HsMetricsCalculator instance to do the real work.
     */
    protected int doWork() {
        IoUtil.assertFileIsReadable(BAIT_INTERVALS);
        IoUtil.assertFileIsReadable(TARGET_INTERVALS);
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(METRICS_FILE);

        HsMetricsCalculator calculator = new HsMetricsCalculator(BAIT_INTERVALS, TARGET_INTERVALS);
        SAMFileReader sam = new SAMFileReader(INPUT);
        calculator.analyze(sam.iterator());

        MetricsFile<HsMetrics, Integer> metrics = getMetricsFile();
        metrics.addMetric(calculator.getMetrics());

        metrics.write(METRICS_FILE);
        return 0;
    }
}
