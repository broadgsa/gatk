package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A parallelizable walker designed to quickly aggregate relevant coverage statistics across samples in the input
 * file. Assesses the mean and median granular coverages of each sample, and generates part of a cumulative
 * distribution of % bases and % targets covered for certain depths. The granularity of DOC can be set by command
 * line arguments.
 *
 * // todo -- allow for user to set linear binning (default is logarithmic)
 * // todo -- add per target (e.g. regional) aggregation
 *
 * @Author chartl
 * @Date Feb 22, 2010
 */
@By(DataSource.REFERENCE)
public class CoverageStatistics extends LocusWalker<Map<String,Integer>, DepthOfCoverageStats> implements TreeReducible<DepthOfCoverageStats> {
    @Argument(fullName = "start", doc = "Starting (left endpoint) for granular binning", required = false)
    int start = 1;
    @Argument(fullName = "stop", doc = "Ending (right endpoint) for granular binning", required = false)
    int stop = 1000;
    @Argument(fullName = "nBins", doc = "Number of bins to use for granular binning", required = false)
    int nBins = 20;
    @Argument(fullName = "minMappingQuality", shortName = "mmq", doc = "Minimum mapping quality of reads to count towards depth. Defaults to 50.", required = false)
    byte minMappingQuality = 50;
    @Argument(fullName = "minBaseQuality", shortName = "mbq", doc = "Minimum quality of bases to count towards depth. Defaults to 20.", required = false)
    byte minBaseQuality = 20;
    @Argument(fullName = "perLocusStatisticsFile", shortName = "locusFile", doc = "File to output per-locus statistics to; if unprovided these will not be calculated", required = false)
    File perLocusStatisticsFile = null;
    @Argument(fullName = "perSampleStatisticsFile", shortName = "sampleFile", doc = "File to output per-sample statistics to; if unprovided will go to standard (-o) output", required = false)
    File perSampleStatisticsFile = null;
    @Argument(fullName = "summaryStatisticsFile", shortName = "summaryFile", doc = "File to output summary (mean, median) statistics to; if unprovided will go to standard (-o) output", required = false)
    File summaryStatisticsFile = null;

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD WALKER METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public DepthOfCoverageStats reduceInit() {
        List<Set<String>> samplesByReaders = getToolkit().getSamplesByReaders();
        DepthOfCoverageStats stats = new DepthOfCoverageStats(DepthOfCoverageStats.calculateBinEndpoints(start,stop,nBins));

        for ( Set<String> sampleSet : samplesByReaders ) {
            for ( String sample : sampleSet ) {
                stats.addSample(sample);
            }
        }

        if ( perLocusStatisticsFile != null ) {
            stats.initializeLocusCounts();
        }

        return stats;
    }

    public Map<String,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Map<String,StratifiedAlignmentContext> contextsBySample =
                                            StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
        HashMap<String,Integer> depthBySample = new HashMap<String,Integer>();

        for ( String sample : contextsBySample.keySet() ) {
            AlignmentContext sampleContext = contextsBySample.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);
            int properDepth = 0;
            for ( PileupElement e : sampleContext.getBasePileup() ) {
                if ( e.getQual() >= minBaseQuality && e.getMappingQual() >= minMappingQuality ) {
                    properDepth++;
                }
            }

            depthBySample.put(sample,properDepth);
        }

        return depthBySample;
    }

    public DepthOfCoverageStats reduce(Map<String,Integer> thisMap, DepthOfCoverageStats prevReduce) {
        prevReduce.update(thisMap);
        return prevReduce;
    }

    public DepthOfCoverageStats treeReduce(DepthOfCoverageStats left, DepthOfCoverageStats right) {
        left.merge(right);
        return left;
    }

    public void onTraversalDone(DepthOfCoverageStats coverageProfiles) {
        printSummary(out,summaryStatisticsFile,coverageProfiles);
        printPerSample(out,perSampleStatisticsFile,coverageProfiles);
        printPerLocus(perLocusStatisticsFile,coverageProfiles);
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // HELPER OUTPUT METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    private void printPerSample(PrintStream out, File optionalFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(out,optionalFile);
        int[] leftEnds = stats.getEndpoints();
        StringBuilder hBuilder = new StringBuilder();
        hBuilder.append("\t");
        hBuilder.append(String.format("[0,%d)\t",leftEnds[0]));
        for ( int i = 1; i < leftEnds.length; i++ )
            hBuilder.append(String.format("[%d,%d)\t",leftEnds[i-1],leftEnds[i]));
        hBuilder.append(String.format("[%d,inf)%n",leftEnds[leftEnds.length-1]));
        output.print(hBuilder.toString());
        Map<String,int[]> histograms = stats.getHistograms();
        for ( String s : histograms.keySet() ) {
            StringBuilder sBuilder = new StringBuilder();
            sBuilder.append(String.format("%s",s));
            for ( int count : histograms.get(s) ) {
                sBuilder.append(String.format("\t%d",count));
            }
            sBuilder.append(String.format("%n"));
            output.print(sBuilder.toString());
        }
    }

    private void printPerLocus(File locusFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(null,locusFile);
        if ( output == null ) {
            return;
        }

        int[] endpoints = stats.getEndpoints();
        int samples = stats.getHistograms().size();

        int[][] baseCoverageCumDist = stats.getLocusCounts();

        // rows - # of samples
        // columns - depth of coverage

        StringBuilder header = new StringBuilder();
        for ( int d : endpoints ) {
            header.append(String.format("\t%d",d));
        }
        header.append(String.format("%n"));

        output.print(header);

        for ( int row = samples; row > 0; row ++ ) {
            output.printf("%s_%d\t","NSamples",row);
            for ( int depthBin = 0; depthBin < baseCoverageCumDist[0].length; depthBin ++ ) {
                output.printf("%d\t",baseCoverageCumDist[row][depthBin]);
            }
            output.printf("%n");
        }
    }

    private PrintStream getCorrectStream(PrintStream out, File optionalFile) {
        PrintStream output;
        if ( optionalFile == null ) {
            output = out;
        } else {
            try {
                output = new PrintStream(optionalFile);
            } catch ( IOException e ) {
                logger.warn("Error opening the output file "+optionalFile.getAbsolutePath()+". Defaulting to stdout");
                output = out;
            }
        }
        
        return output;
    }

    private void printSummary(PrintStream out, File optionalFile, DepthOfCoverageStats stats) {
        PrintStream output = getCorrectStream(out,optionalFile);

        output.printf("%s\t%s\t%s\t%s\t%s%n","sample_id","mean","granular_third_quartile","granular_median","granular_first_quartile");

        Map<String,int[]> histograms = stats.getHistograms();
        Map<String,Double> means = stats.getMeans();
        int[] leftEnds = stats.getEndpoints();

        for ( String s : histograms.keySet() ) {
            int median = getQuantile(histograms.get(s),0.5);
            int q1 = getQuantile(histograms.get(s),0.25);
            int q3 = getQuantile(histograms.get(s),0.75);
            output.printf("%s\t%.2f\t%d\t%d\t%d%n",s,means.get(s),leftEnds[q3],leftEnds[median],leftEnds[q1]);
        }

        output.printf("%s\t%.2f\t%s\t%s\t%s%n","Total",means.get(DepthOfCoverageStats.ALL_SAMPLES),"N/A","N/A","N/A");
    }

    private int getQuantile(int[] histogram, double prop) {
        int total = 0;

        for ( int i = 0; i < histogram.length; i ++ ) {
            total += histogram[i];
        }

        int counts = 0;
        int bin = -1;
        while ( counts < prop*total ) {
            counts += histogram[bin+1];
            bin++;
        }

        return bin;
    }
}

class DepthOfCoverageStats {
    public static String ALL_SAMPLES = "ALL_COMBINED_SAMPLES";
    // STANDARD (constantly updated) DATA
    private Map<String,int[]> granularHistogramBySample; // holds the counts per each bin
    private Map<String,Double> meanCoverages; // holds mean coverage per sample
    private int[] binLeftEndpoints; // describes the left endpoint for each bin
    private int[][] locusCoverageCounts; // holds counts of number of bases with >=X samples at >=Y coverage
    private boolean tabulateLocusCounts = false;
    private int nLoci; // number of loci seen
    // TEMPORARY DATA (not worth re-instantiating every time)
    private int[] locusHistogram; // holds a histogram for each locus; reset after each update() call
    private int totalDepth; // holds the total depth of coverage for each locus; reset after each update() call

    ////////////////////////////////////////////////////////////////////////////////////
    // STATIC METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public static int[] calculateBinEndpoints(int lower, int upper, int bins) {
        if ( bins > upper - lower || lower < 1 ) {
            throw new IllegalArgumentException("Illegal argument to calculateBinEndpoints; "+
                    "lower bound must be at least 1, and number of bins may not exceed stop - start");
        }

        int[] binLeftEndpoints = new int[bins+1];
        binLeftEndpoints[0] = lower;

        int length = upper - lower;
        double scale = Math.log10((double) length)/bins;

        for ( int b = 1; b < bins ; b++ ) {
            int leftEnd = lower + (int) Math.floor(Math.pow(10.0,(b-1.0)*scale));

            while ( leftEnd <= binLeftEndpoints[b-1] ) {
                leftEnd++;
            }

            binLeftEndpoints[b] = leftEnd;
        }

        binLeftEndpoints[binLeftEndpoints.length-1] = upper;

        return binLeftEndpoints;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // INITIALIZATION METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public DepthOfCoverageStats(int[] leftEndpoints) {
        this.binLeftEndpoints = leftEndpoints;
        granularHistogramBySample = new HashMap<String,int[]>();
        meanCoverages = new HashMap<String,Double>();
        meanCoverages.put(DepthOfCoverageStats.ALL_SAMPLES,0.0);
        nLoci = 0;
        totalDepth = 0;
    }

    public void addSample(String sample) {
        if ( granularHistogramBySample.containsKey(sample) ) {
            return;
        }

        int[] binCounts = new int[this.binLeftEndpoints.length];
        for ( int b = 0; b < binCounts.length; b ++ ) {
            binCounts[b] = 0;
        }

        granularHistogramBySample.put(sample,binCounts);
        meanCoverages.put(sample,0.0);
    }

    public void initializeLocusCounts() {
        locusCoverageCounts = new int[granularHistogramBySample.size()][binLeftEndpoints.length];
        locusHistogram = new int[binLeftEndpoints.length];
        for ( int b = 0; b < binLeftEndpoints.length; b ++ ) {
            for ( int a = 0; a < granularHistogramBySample.size(); a ++ ) {
                locusCoverageCounts[a][b] = 0;
            }
            locusHistogram[b] = 0;
        }

        tabulateLocusCounts = true;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // UPDATE METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public void update(Map<String,Integer> depthBySample) {
        int b;
        for ( String sample : granularHistogramBySample.keySet() ) {
            if ( depthBySample.containsKey(sample) ) {
                b = updateSample(sample,depthBySample.get(sample));
                totalDepth += depthBySample.get(sample);
            } else {
                b = updateSample(sample,0);
            }

            if ( tabulateLocusCounts ) {
                for ( int i = 0; i <= b; i ++ ) {
                    locusHistogram[i]++;
                }
            }
        }

        double meanDepth = meanCoverages.get(DepthOfCoverageStats.ALL_SAMPLES);
        double newMean = ( meanDepth*nLoci + (double) totalDepth )/( nLoci + 1 );
        meanCoverages.put(DepthOfCoverageStats.ALL_SAMPLES,newMean);
        updateLocusCounts(locusHistogram);

        nLoci++;
        totalDepth = 0;
    }

    private int updateSample(String sample, int depth) {
        double mean = meanCoverages.get(sample);
        double newMean = ( nLoci*mean + (double) depth )/(nLoci + 1.0);
        meanCoverages.put(sample,newMean);

        int[] granularBins = granularHistogramBySample.get(sample);
        for ( int b = 1; b < granularBins.length; b ++ ) {
            if ( depth < binLeftEndpoints[b] ) {
                granularBins[b-1]++;
                return b ;
            }
        }

        granularBins[granularBins.length-1]++; // greater than all left-endpoints
        return granularBins.length-1;
    }

    public void merge(DepthOfCoverageStats newStats) {
        this.mergeSamples(newStats);
        if ( this.tabulateLocusCounts && newStats.tabulateLocusCounts ) {
            this.mergeLocusCounts(newStats.getLocusCounts());
        }

        double totalMean = (meanCoverages.get(DepthOfCoverageStats.ALL_SAMPLES)*nLoci +
                newStats.getMeans().get(DepthOfCoverageStats.ALL_SAMPLES)*newStats.getTotalLoci()) /
                ( nLoci + newStats.getTotalLoci());

        meanCoverages.put(DepthOfCoverageStats.ALL_SAMPLES,totalMean);
        nLoci += newStats.getTotalLoci();
    }

    private void mergeSamples(DepthOfCoverageStats otherStats) {
        Map<String,int[]> otherHistogram = otherStats.getHistograms();
        Map<String,Double> otherMeans = otherStats.getMeans();
        for ( String s : granularHistogramBySample.keySet() ) {
            int[] internalCounts = granularHistogramBySample.get(s);
            int[] externalCounts = otherHistogram.get(s);
            for ( int b = 0; b < internalCounts.length; b++ ) {
                internalCounts[b] += externalCounts[b];
            }

            double internalMean = meanCoverages.get(s);
            double externalMean = otherMeans.get(s);
            double newMean = ( internalMean*nLoci + externalMean*otherStats.getTotalLoci())/(nLoci+otherStats.getTotalLoci());

            meanCoverages.put(s,newMean);
        }
    }

    private void mergeLocusCounts( int[][] otherCounts ) {
        for ( int a = 0; a < locusCoverageCounts.length; a ++ ) {
            for ( int b = 0; b < locusCoverageCounts[0].length; b ++ ) {
                locusCoverageCounts[a][b] += otherCounts[a][b];
            }
        }
    }

    /*
     * Update locus counts -- takes an array in which the number of samples
     * with depth ABOVE [i] is held. So if the bin left endpoints were 2, 5, 10
     * then we'd have an array that represented:
     * [# samples with depth 0 - inf], [# samples with depth 2 - inf],
     * [# samples with depth 5 - inf], [# samples with depth 10-inf];
     *
     * this is
     * @argument cumulativeSamplesByDepthBin - see above
     */
    private void updateLocusCounts(int[] cumulativeSamplesByDepthBin) {
        if ( tabulateLocusCounts ) {
            for ( int bin = 0; bin < cumulativeSamplesByDepthBin.length; bin ++ ) {
                int numSamples = cumulativeSamplesByDepthBin[bin];
                for ( int i = 0; i < numSamples; i ++ ) {
                    locusCoverageCounts[i][bin]++;
                }

                cumulativeSamplesByDepthBin[bin] = 0; // reset counts in advance of next update()
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // ACCESSOR METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public Map<String,int[]> getHistograms() {
        return granularHistogramBySample;
    }

    public int[][] getLocusCounts() {
        return locusCoverageCounts;
    }

    public int[] getEndpoints() {
        return binLeftEndpoints;
    }

    public Map<String,Double> getMeans() {
        return meanCoverages;
    }

    public int getTotalLoci() {
        return nLoci;
    }

}
