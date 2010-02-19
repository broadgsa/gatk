package org.broadinstitute.sting.oneoffprojects.firehosesummary;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.*;
import java.util.*;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 11, 2010
 */

class FirehoseSummaryCLP extends CommandLineProgram {
    @Argument(fullName = "depthOfCoverageFile", shortName = "doc", doc="Path to the depth of coverage file", required=true)
    private File depthOfCoverage = null;
//    @Argument(fullName = "contaminationFile", shortName = "con", doc="Path to the contamination file", required=true)
//    private File contamination = null;
//    @Argument(fullName = "errorRateFile", shortName = "err", doc="Path to the error rate file", required=true)
//    private File errorRate = null;
//    @Argument(fullName = "zipFiles", shortName = "zip", doc="List of paths to zip files which contain summary metrics files", required=false)
//    private String zipFiles = null;

    private static String R_SCRIPT = "plotFirehoseDataQCMetrics.R";
    private static String SCRIPT_DOC_FLAG = "DOC";

    protected int execute() {
//        SummaryFileCollection metricsFiles = getFileHandles();
        List<DepthStatisticsCalculator> depthStats = calculateDepthStatistics(depthOfCoverage);
        String docSummary = makeDOCPlots(depthStats);
        return 1;
    }

    private String makeDOCPlots(List<DepthStatisticsCalculator> calcs) {
        StringBuilder summaryBuilder = new StringBuilder();
        int numSamplesWithGoodBaseCoverage=0;
        int numSamplesWithGoodMeanCoverage = 0;
        int numSamplesWithGoodTargetCoverage=0;
        double aggregateMeanCoverage = -1d;
        double aggregateStdDevCoverage = -1d;
        double aggregateSkewCoverage = -1d;
        PrintWriter locusOutput;
        PrintWriter targetOutput;
        File locusOutputFile;
        File targetOutputFile;
        try {
            locusOutputFile = File.createTempFile("locus_output","txt");
            locusOutput = new PrintWriter(locusOutputFile);
            targetOutputFile = File.createTempFile("target_output","txt");
            targetOutput = new PrintWriter(targetOutputFile);
        } catch ( IOException e ) {
            throw new StingException("Error creating temporary output files for R-plotting. Perhaps user permissions do not allow creation of files.",e);
        }

        for ( DepthStatisticsCalculator calc : calcs ) {

            if ( calc.getName().equalsIgnoreCase("total_coverage") || calc.getName().equalsIgnoreCase("coverage_without_deletions") ) {
                // update summary statistics
                if ( calc.getName().equalsIgnoreCase("total_coverage") ) {
                    aggregateMeanCoverage = calc.getMean();
                    aggregateStdDevCoverage = Math.sqrt(calc.getVar());
                    aggregateSkewCoverage = calc.getSkew();
                }
            } else {
                // update sample-based summary statistics and print output to file
                locusOutput.printf("%s\t%f\t%f\t%f\t%f\t%f\t%f%n",calc.getName(),calc.getLocusProportions());
                targetOutput.printf("%s\t%f\t%f\t%f\t%f\t%f\t%f%n", calc.getName(),calc.getTargetProportions());
                if ( calc.getPercentWellCoveredLoci() > 80 ) {
                    numSamplesWithGoodBaseCoverage++;
                }
                if ( calc.getPercentWellCoveredTargets() > 90 ) {
                    numSamplesWithGoodTargetCoverage++;
                }
                if ( calc.getMean() > 10 ) {
                    numSamplesWithGoodMeanCoverage++;
                }
            }
        }

        //invokeRScript(R_SCRIPT,SCRIPT_DOC_FLAG,"PercentOfBasesCoveredBySample",locusOutputFile.getAbsolutePath());
        //invokeRScript(R_SCRIPT,SCRIPT_DOC_FLAG,"PercentOfTargetsCoveredBySample",targetOutputFile.getAbsolutePath());
        return "temporary";
    }

//    private SummaryFileCollection getFileHandles() {
//        if ( zipFiles == null ) {
//            return null;
//        }
//
//        SummaryFileCollection summaryFiles = new SummaryFileCollection();
//        for ( String zipFile : zipFiles.split(",") ) {
//            summaryFiles.process(zipFile);
//        }
//
//        return summaryFiles;
//    }

    private List<DepthStatisticsCalculator> calculateDepthStatistics(File docFile) {
        BufferedReader docReader;

        try {
            docReader = new BufferedReader( new FileReader(docFile) );
        } catch ( IOException e) {
            throw new StingException("The file "+docFile.getAbsolutePath()+" could not be opened...",e);
        }

        String locusHeader = getDOCSectionHeader(docReader);
        List<DepthStatisticsCalculator> docCalculators = instantiateDOCCalculators(locusHeader);

        updateLocusInfo(docCalculators,docReader);

        String targetHeader = getDOCSectionHeader(docReader);

        updateTargetInfo(docCalculators,docReader);

        return docCalculators;
    }

    private void updateLocusInfo(List<DepthStatisticsCalculator> calcs, BufferedReader reader) {

        String docLocLine;
        try {
            docLocLine = reader.readLine();
            while ( ! isEndOfSection(docLocLine) ) {
                int offset = -1;
                for ( String entry : docLocLine.split("\t") ) {
                    if ( offset > -1 ) {
                        calcs.get(offset).updateLocus(Integer.parseInt(entry));
                    }
                    offset++;
                }
            }
        } catch ( IOException e) {
            throw new StingException("Error reading locus depth of coverage information",e);
        }

    }

    private void updateTargetInfo(List<DepthStatisticsCalculator> calcs, BufferedReader reader) {

        String docLocLine;
        try {
            docLocLine = reader.readLine();
            while ( ! isEndOfSection(docLocLine) ) {
                int offset = -1;
                int targetSize = 0;
                for ( String entry : docLocLine.split("\t") ) {
                    if ( offset == -1 ) {
                        targetSize = parseInterval(entry);
                    } else {
                        calcs.get(offset).updateTargets(targetSize,Integer.parseInt(entry));
                    }
                    offset++;
                }
            }
        } catch ( IOException e ) {
            throw new StingException("Error reading target depth of coverage information",e);
        }

    }

    private int parseInterval(String interval) {
        String startstop = interval.split(":")[1];
        int start = Integer.parseInt(startstop.split("-")[0]);
        int stop = Integer.parseInt(startstop.split("-")[1]);
        return stop - start;
    }

    private boolean isDOCSectionSeparator( String line ) {
        return line.contains("_COVERAGE_SECTION");
    }

    private boolean isEndOfSection( String line ) {
        // sections delimited by empty line
        return line.equalsIgnoreCase("");
    }

    private String getDOCSectionHeader(BufferedReader reader) {
        String header;
        try {
            do {
                header = reader.readLine();
            } while ( ! isDOCSectionSeparator(header) );

            header = reader.readLine();

        } catch (IOException e) {
            throw new StingException("Error reading depth of coverage file",e);
        }

        return header;
    }

    private List<DepthStatisticsCalculator> instantiateDOCCalculators(String header) {
        List<DepthStatisticsCalculator> calcs = new ArrayList<DepthStatisticsCalculator>();

        int offset = -1;
        for ( String entry : header.split("\t") ) {
            if ( offset > -1 ) {
                calcs.add(new DepthStatisticsCalculator(entry));
            }
            offset++;
        }

        return calcs;
    }

}

public class GenerateFirehoseSummary {

    public static void main(String[] args) {
        FirehoseSummaryCLP clp = new FirehoseSummaryCLP();
        CommandLineProgram.start(clp, args);
        System.exit(0);
    }
}

//class SummaryFileCollection {
//
//    // container class for files we'll be summarizing
//
//    public Map<String,File> fingerprintSummaryFiles;
//    public Map<String,File> hybridSelectionMetricsFiles;
//    public Map<String,File> insertSizeDistributionFiles;
//    public Map<String,File> alignmentMetricsFiles;
//
//    public SummaryFileCollection() {
//        fingerprintSummaryFiles = new HashMap<String,File>();
//        hybridSelectionMetricsFiles = new HashMap<String, File>();
//        insertSizeDistributionFiles = new HashMap<String,File>();
//        alignmentMetricsFiles = new HashMap<String,File>();
//    }
//
//    public void process(String zipFilePath) {
//        String sampleName = zipFilePath.split("_sequencing_metrics.zip")[0].split("_")[1];
//        File fingerprintSummaryFile = new File(sampleName+".summary_fingerprint_metrics");
//        File hybridSelectionFile = new File(sampleName+".hybrid_selection_metrics");
//        File insertSizeFile = new File(sampleName+".insert_size_metrics");
//        File alignmentFile = new File(sampleName+".alignment_metrics");
//
//        String command = "unzip "+zipFilePath;
//        try {
//            Process p = Runtime.getRuntime().exec(command);
//        } catch (IOException e) {
//            throw new RuntimeException("Could not unzip the file "+zipFilePath);
//        }
//
//        fingerprintSummaryFiles.put(sampleName,fingerprintSummaryFile);
//        hybridSelectionMetricsFiles.put(sampleName,hybridSelectionFile);
//        insertSizeDistributionFiles.put(sampleName,insertSizeFile);
//        alignmentMetricsFiles.put(sampleName,alignmentFile);
//    }
//}