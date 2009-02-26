/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.sam;

import java.io.File;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.io.IoUtil;
import edu.mit.broad.picard.metrics.MetricsFile;
import edu.mit.broad.picard.util.Histogram;
import edu.mit.broad.picard.util.Log;
import edu.mit.broad.picard.util.RExecutor;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMRecord;
import edu.mit.broad.sam.util.CloseableIterator;

/**
 * Command line program to read non-duplicate insert sizes, create a histogram
 * and report distribution statistics.
 *
 * @author Doug Voet
 */
public class CollectInsertSizeMetrics extends CommandLineProgram {
    private static Log log = Log.getInstance(CollectInsertSizeMetrics.class);
    private static final String HISTOGRAM_R_SCRIPT = "edu/mit/broad/picard/sam/insertSizeHistogram.R";
    // Usage and parameters
    @Usage(programVersion="1.0") 
    public String USAGE = "Reads a SAM or BAM file and writes a file containing metrics about " +
    		"the statistical distribution of insert size (excluding duplicates) " +
    		"and generates a histogram plot.\n";
    @Option(shortName="I", doc="SAM or BAM file") public File INPUT;
    @Option(shortName="O", doc="File to write insert size metrics to") public File OUTPUT;
    @Option(shortName="H", doc="File to write insert size histogram chart to") public File HISTOGRAM_FILE;
	
    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new CollectInsertSizeMetrics().instanceMain(argv));
    }

	@Override
	protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(HISTOGRAM_FILE);

        SAMFileReader in = new SAMFileReader(INPUT);
        MetricsFile<InsertSizeMetrics, Integer> file = collectMetrics(in.iterator());
		in.close();
        
        file.write(OUTPUT);
        
        if (file.getMetrics().get(0).READ_PAIRS == 0) {
            log.warn("Input file did not contain any records with insert size information.");
        } else  {
            int rResult = RExecutor.executeFromClasspath(
                    HISTOGRAM_R_SCRIPT, 
                    OUTPUT.getAbsolutePath(), 
                    HISTOGRAM_FILE.getAbsolutePath(),
                    INPUT.getName());

            if (rResult != 0) {
                throw new PicardException("R script " + HISTOGRAM_R_SCRIPT + " failed with return code " + rResult);
            }
        }
        
        return 0;
	}

	/**
	 * Does all the work of iterating through the sam file and collecting insert size metrics.
	 */
	MetricsFile<InsertSizeMetrics, Integer> collectMetrics(CloseableIterator<SAMRecord> samIterator) {
		Histogram<Integer> insertSizeHistogram = new Histogram<Integer>("insert_size", "count");
        while (samIterator.hasNext()) {
			SAMRecord record = samIterator.next();
			if (skipRecord(record)) {
				continue;
			}
			
			int insertSize = Math.abs(record.getInferredInsertSize());
            insertSizeHistogram.increment(insertSize);
		}

        MetricsFile<InsertSizeMetrics, Integer> file = new MetricsFile<InsertSizeMetrics, Integer>();
        file.setHistogram(insertSizeHistogram);
        InsertSizeMetrics metrics = new InsertSizeMetrics();
        metrics.READ_PAIRS = (long) insertSizeHistogram.getCount();
        metrics.MAX_INSERT_SIZE = (int) insertSizeHistogram.getMax();
        metrics.MIN_INSERT_SIZE = (int) insertSizeHistogram.getMin();
        metrics.MEAN_INSERT_SIZE = insertSizeHistogram.getMean();
        metrics.STANDARD_DEVIATION = insertSizeHistogram.getStandardDeviation();
        metrics.MEDIAN_INSERT_SIZE = insertSizeHistogram.getMedian();

        final double total   = insertSizeHistogram.getCount();
        final double median  = insertSizeHistogram.getMedian();
        double covered = 0;
        double low  = median;
        double high = median;

        while (low >= insertSizeHistogram.getMin() || high <= insertSizeHistogram.getMax()) {
            Histogram<Integer>.Bin lowBin = insertSizeHistogram.get((int) low);
            if (lowBin != null) covered += lowBin.getValue();

            if (low != high) {
                Histogram<Integer>.Bin highBin = insertSizeHistogram.get((int) high);
                if (highBin != null) covered += highBin.getValue();
            }

            double percentCovered = covered / total;
            int distance = (int) (high - low) + 1;
            if (percentCovered >= 0.1  && metrics.WIDTH_OF_10_PERCENT == 0) metrics.WIDTH_OF_10_PERCENT = distance;
            if (percentCovered >= 0.2  && metrics.WIDTH_OF_20_PERCENT == 0) metrics.WIDTH_OF_20_PERCENT = distance;
            if (percentCovered >= 0.3  && metrics.WIDTH_OF_30_PERCENT == 0) metrics.WIDTH_OF_30_PERCENT = distance;
            if (percentCovered >= 0.4  && metrics.WIDTH_OF_40_PERCENT == 0) metrics.WIDTH_OF_40_PERCENT = distance;
            if (percentCovered >= 0.5  && metrics.WIDTH_OF_50_PERCENT == 0) metrics.WIDTH_OF_50_PERCENT = distance;
            if (percentCovered >= 0.6  && metrics.WIDTH_OF_60_PERCENT == 0) metrics.WIDTH_OF_60_PERCENT = distance;
            if (percentCovered >= 0.7  && metrics.WIDTH_OF_70_PERCENT == 0) metrics.WIDTH_OF_70_PERCENT = distance;
            if (percentCovered >= 0.8  && metrics.WIDTH_OF_80_PERCENT == 0) metrics.WIDTH_OF_80_PERCENT = distance;
            if (percentCovered >= 0.9  && metrics.WIDTH_OF_90_PERCENT == 0) metrics.WIDTH_OF_90_PERCENT = distance;
            if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0) metrics.WIDTH_OF_99_PERCENT = distance;

            --low;
            ++high;
        }

        file.addMetric(metrics);
        
		return file;
	}

    /**
     * Figures out whether or not the record should be included in the counting of insert sizes
     */
    private boolean skipRecord(SAMRecord record) {
        return !record.getReadPairedFlag() || 
                record.getMateUnmappedFlag() || 
                record.getFirstOfPairFlag() || 
                record.getNotPrimaryAlignmentFlag() || 
                record.getDuplicateReadFlag() ||
                record.getInferredInsertSize() == 0;
    }

}
