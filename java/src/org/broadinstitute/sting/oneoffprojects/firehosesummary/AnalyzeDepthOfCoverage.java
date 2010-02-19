package org.broadinstitute.sting.oneoffprojects.firehosesummary;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.cmdLine.CommandLineProgram;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 18, 2010
 */

class AnalyzeDepthCLP extends CommandLineProgram {
    @Argument(fullName = "depthOfCoverageFile", shortName = "df", doc = "The Depth of Coverage output file", required = true)
    public File docFile = null;
    @Argument(fullName = "summaryFile", shortName = "sf", doc = "The summary file to which to output", required = true)
    public File summaryFile = null;
    @Argument(fullName = "plotBaseName", shortName = "bn", doc = "The base name for the plot files (e.g. 'foo' yields plots 'foo_DoC_by_sample.pdf'). Please ensure this name contains no spaces.", required = false)
    public String plotBaseName = "DepthAnalysis";
    @Argument(fullName = "pathToRScript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "./";

    private boolean containsByLocus = false;
    private boolean containsByTarget = false;

    ///////////////////////////////////////////////////////////////////////////////////
    // CONSTANT VALUES: SUMMARY STRING FOR NO INFORMATION, R-SCRIPT ARGUMENTS, ETC
    ///////////////////////////////////////////////////////////////////////////////////

    private final String DEFAULT_SUMMARY_STRING = "No Summary Information";
    private final String PER_LOCUS_R_ARGUMENTS = "PlotInterleavedRows depth_of_coverage\\;proportion_of_bases_above\\;Per_Sample_Depth_of_Coverage\\;"+plotBaseName+"_per_locus";
    private final String PER_TARGET_R_ARGUMENTS = "PlotInterleavedRows depth_of_coverage\\;proportion_of_targets_with_mean_coverage_above\\;Per_Sample_Average_DoC_Over_Targets\\;"+plotBaseName+"_per_target";

    ///////////////////////////////////////////////////////////////////////////////////
    // ANALYSIS START: CALCULATE STATISTICS, WRITE IN R-READABLE FORMAT, MAKE PLOTS
    ///////////////////////////////////////////////////////////////////////////////////

    protected int execute() {
        List<DepthStatisticsCalculator> depthStats = calculateDepthStatistics(docFile);
        String perLocusSummary = DEFAULT_SUMMARY_STRING;
        String perTargetSummary = DEFAULT_SUMMARY_STRING;

        if ( containsByLocus ) {
            File baseSummaryTable = writeBaseSummaryFile(depthStats);
            perLocusSummary = generatePerLocusSummary(baseSummaryTable,depthStats);
        }

        if ( containsByTarget ) {
            File targetSummaryTable = writeTargetSumamryFile(depthStats);
            perTargetSummary = generatePerTargetSummary(targetSummaryTable, depthStats);
        }

        writeSummaryInfoFile(summaryFile,perLocusSummary,perTargetSummary);

        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // OPEN AND WRITE FINAL SUMMARY DOC FILE
    ///////////////////////////////////////////////////////////////////////////////////

    private void writeSummaryInfoFile(File sFile, String locusSummary, String targetSummary) {
        PrintWriter writer;
        try {
            writer = new PrintWriter(sFile);
            writer.printf("%s%n","##Depth of coverage summary file");
            writer.printf("%s%n","##Well_Covered_Samples_By_Base - % of samples with >80% bases covered to 10x");
            writer.printf("%s%n","##Well_Covered_Samples_By_Mean - % of samples with mean coverage > 10x");
            writer.printf("%s%n%n","##Well_Covered_Samples_By_Target - % of samples with >80% targets covered to 10x");
            if ( containsByLocus )
                writer.printf("%s%n",locusSummary);
            if ( containsByTarget )
                writer.printf("%s",targetSummary);
            writer.close();
        } catch (IOException e) {
            throw new StingException("Error writing final depth of coverage summary file",e);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // CALL R-SCRIPTS AND GENERATE OVERALL SUMMARY FILES
    ///////////////////////////////////////////////////////////////////////////////////

    private String generatePerLocusSummary(File rReadablePlotFile, List<DepthStatisticsCalculator> calcs) {
        String rCommand = PATH_TO_RSCRIPT+" "+PATH_TO_RESOURCES+" "+rReadablePlotFile.getAbsolutePath()+" "+PER_LOCUS_R_ARGUMENTS;
        try {
            Process p = Runtime.getRuntime().exec(rCommand);
        } catch ( IOException e ) {
            throw new StingException("Error executing r command for per locus plot generation",e);
        }

        StringBuilder summary = new StringBuilder();
        summary.append(String.format("%s%n","PER_LOCUS_SUMMARY"));
        int numSamples = calcs.size()-2;
        int numGoodSamples = 0;
        int numGoodSamplesByMeanCvg = 0;
        double totalAvgCoverage = -1;
        double totalStdevCoverage = -1;

        for ( DepthStatisticsCalculator calc : calcs ) {
            if ( calc.getName().equalsIgnoreCase("total_coverage")) {
                totalAvgCoverage = calc.getMean();
                totalStdevCoverage = Math.sqrt(calc.getVar());
            } else if ( ! calc.getName().equalsIgnoreCase("coverage_without_deletions") ) {
                if ( calc.getPercentWellCoveredLoci() > 0.8 ) {
                    numGoodSamples++;
                }

                if ( calc.getMean() > 10 ) {
                    numGoodSamplesByMeanCvg++;
                }
            }
        }

        summary.append(String.format("%s\t%f%n","Average_Coverage:",totalAvgCoverage));
        summary.append(String.format("%s\t%f%n","Stdev_Coverage:",totalStdevCoverage));
        summary.append(String.format("%s\t%.2f%n","%Well_Covered_Samples_By_Base", ( (double) numGoodSamples*100 )/( (double) numSamples)));
        summary.append(String.format("%s\t%.2f%n","%Well_Covered_Samples_By_Mean", ( (double) numGoodSamplesByMeanCvg*100) / ( (double) numSamples )));

        return summary.toString();
    }

    private String generatePerTargetSummary(File rReadablePlotFile, List<DepthStatisticsCalculator> calcs) {
        String rCommand = PATH_TO_RSCRIPT+" "+PATH_TO_RESOURCES+" "+rReadablePlotFile.getAbsolutePath()+" "+PER_TARGET_R_ARGUMENTS;
        try {
            Process p = Runtime.getRuntime().exec(rCommand);
        } catch ( IOException e ) {
            throw new StingException("Error executing r command for per locus plot generation",e);
        }

        StringBuilder summary = new StringBuilder();
        summary.append(String.format("%s%n","PER_TARGET_SUMMARY"));
        int numSamples = calcs.size()-2;
        int numGoodSamples = 0;

        for ( DepthStatisticsCalculator calc : calcs ) {
            if ( calc.getName().equalsIgnoreCase("total_coverage")) {
                // do nothing
            } else if ( ! calc.getName().equalsIgnoreCase("coverage_without_deletions") ) {
                if ( calc.getPercentWellCoveredTargets() > 0.8 ) {
                    numGoodSamples++;
                }
            }
        }

        summary.append(String.format("%s\t%.2f%n","%Well_Covered_Samples_By_Target", ( (double) numGoodSamples*100) / ( (double) numSamples )));

        return summary.toString();
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // R-READABLE TEMPORARY FILE CREATION
    ///////////////////////////////////////////////////////////////////////////////////

    private File writeBaseSummaryFile(List<DepthStatisticsCalculator> calcs) {
        File perLocusSummaryFile;

        try {
            perLocusSummaryFile = File.createTempFile(plotBaseName+"_per_locus_summary",".txt");
        } catch ( IOException e ) {
            throw new StingException("Could not create a temporary file. Please check the permissions of the directory you are running in, and that the base name is not a filepath.",e);
        }

        PrintWriter locusWriter;

        try {
            locusWriter = new PrintWriter(perLocusSummaryFile);
        } catch ( IOException e ) {
            throw new StingException("Locus summary temporary file was created but could not be opened.",e);
        }

        for ( DepthStatisticsCalculator calc : calcs ) {
            if ( ! calc.getName().equalsIgnoreCase("total_coverage") && ! calc.getName().equalsIgnoreCase("coverage_without_deletions") ) {
                locusWriter.printf("%s\t%f\t%f\t%f\t%f\t%f\t%f",calc.getName(),calc.getLocusProportions());
                locusWriter.printf("%s\t%d\t%d\t%d\t%d\t%d\t%d",calc.getName(),calc.getEvalPoints());
            }
        }

        locusWriter.close();
        return perLocusSummaryFile;
    }

    private File writeTargetSumamryFile(List<DepthStatisticsCalculator> calcs) {
        File perTargetSummaryFile;

        try {
            perTargetSummaryFile = File.createTempFile(plotBaseName+"_per_target_summary",".txt");
        } catch ( IOException e ) {
            throw new StingException("Could not create a temporary file. Please check the permissions of the directory you are running in, and that the base name is not a filepath.",e);
        }

        PrintWriter targetWriter;

        try {
            targetWriter = new PrintWriter(perTargetSummaryFile);
        } catch ( IOException e ) {
            throw new StingException("Target summary temporary file was created but could not be opened.",e);
        }

        for ( DepthStatisticsCalculator calc : calcs ) {
            if ( ! calc.getName().equalsIgnoreCase("total_coverage") && ! calc.getName().equalsIgnoreCase("coverage_without_deletions") ) {
                targetWriter.printf("%s\t%f\t%f\t%f\t%f\t%f\t%f",calc.getName(),calc.getTargetProportions());
                targetWriter.printf("%s\t%d\t%d\t%d\t%d\t%d\t%d",calc.getName(),calc.getEvalPoints());
            }
        }

        targetWriter.close();
        return perTargetSummaryFile;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    // READING THE DEPTH OF COVERAGE FILE INTO CALCULATOR OBJECTS
    ///////////////////////////////////////////////////////////////////////////////////

    private List<DepthStatisticsCalculator> calculateDepthStatistics(File docFile) {
        BufferedReader docReader;

        try {
            docReader = new BufferedReader( new FileReader(docFile) );
        } catch ( IOException e) {
            throw new StingException("The file "+docFile.getAbsolutePath()+" could not be opened...",e);
        }

        String locusHeader = getDOCSectionHeader(docReader); // this will read to the first section header
        List<DepthStatisticsCalculator> docCalculators;
        if ( locusHeader != null && locusHeader.equalsIgnoreCase("PER_LOCUS_COVERAGE_SECTION")) {
            containsByLocus = true;
            docCalculators = instantiateDOCCalculators(docReader);
            updateLocusInfo(docCalculators,docReader);
            String targetHeader = getDOCSectionHeader(docReader);
            if ( targetHeader != null && targetHeader.equalsIgnoreCase("PER_TARGET_COVERAGE_SECTION") ) {
                containsByTarget = true;
                updateTargetInfo(docCalculators,docReader);
            } else {
                containsByTarget = false;
            }
        } else if ( locusHeader != null && locusHeader.equalsIgnoreCase("PER_TARGET_COVERAGE_SECTION") ) {
            containsByTarget = true;
            containsByLocus = false;
            docCalculators = instantiateDOCCalculators(docReader);
            updateTargetInfo(docCalculators,docReader);
        } else {
            containsByLocus = false;
            containsByTarget = false;
            docCalculators = null;
        }

        return docCalculators;
    }

    private List<DepthStatisticsCalculator> instantiateDOCCalculators(BufferedReader reader) {
        String header;
        try {
            header = reader.readLine();
        } catch (IOException e) {
            throw new StingException("Unable to read the section header",e);
        }

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

    ///////////////////////////////////////////////////////////////////////////////////
    // FILE IO METHODS -- DEPEND ON DEPTH OF COVERAGE FILE FORMAT
    ///////////////////////////////////////////////////////////////////////////////////

    private boolean isEndOfSection( String line ) {
        // sections delimited by empty line
        return line.equalsIgnoreCase("");
    }

    private String getDOCSectionHeader(BufferedReader reader) {
        String header;
        try {
            do {
                header = reader.readLine();
            } while ( ! isDOCSectionSeparator(header) && header != null);

        } catch (IOException e) {
            throw new StingException("Error reading depth of coverage file",e);
        }

        return header;
    }

    private boolean isDOCSectionSeparator( String line ) {
        return line.contains("_COVERAGE_SECTION");
    }

    private int parseInterval(String interval) {
        String startstop = interval.split(":")[1];
        int start = Integer.parseInt(startstop.split("-")[0]);
        int stop = Integer.parseInt(startstop.split("-")[1]);
        return stop - start;
    }

}

///////////////////////////////////////////////////////////////////////////////////
// PROGRAM START -- THE MAIN() METHOD AND WRAPPER CLASS
///////////////////////////////////////////////////////////////////////////////////

public class AnalyzeDepthOfCoverage {

    public static void main(String[] args) {
        AnalyzeDepthCLP depthAnalysis = new AnalyzeDepthCLP();
        CommandLineProgram.start(depthAnalysis,args);
        System.exit(0);
    }
}
