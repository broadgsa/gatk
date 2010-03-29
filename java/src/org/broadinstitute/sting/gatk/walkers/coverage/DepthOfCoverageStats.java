package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.utils.BaseUtils;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 26, 2010
 */
public class DepthOfCoverageStats {
    ////////////////////////////////////////////////////////////////////////////////////
    // STATIC DATA
    ////////////////////////////////////////////////////////////////////////////////////

    /* none so far */

    ////////////////////////////////////////////////////////////////////////////////////
    // STANDARD DATA
    ////////////////////////////////////////////////////////////////////////////////////

    private Map<String,int[]> granularHistogramBySample; // holds the counts per each bin
    private Map<String,Long> totalCoverages; // holds total coverage per sample
    private int[] binLeftEndpoints; // describes the left endpoint for each bin
    private int[][] locusCoverageCounts; // holds counts of number of bases with >=X samples at >=Y coverage
    private boolean tabulateLocusCounts = false;
    private long nLoci; // number of loci seen
    private long totalDepthOfCoverage;
    private boolean includeDeletions = false;

    ////////////////////////////////////////////////////////////////////////////////////
    // TEMPORARY DATA ( not worth re-instantiating )
    ////////////////////////////////////////////////////////////////////////////////////

    private int[] locusHistogram; // holds a histogram for each locus; reset after each update() call
    private int totalLocusDepth; // holds the total depth of coverage for each locus; reset after each update() call

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
            // todo -- simplify to length^(scale/bins); make non-constant to put bin ends in more "useful"
            // todo -- positions on the number line
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
        totalCoverages = new HashMap<String,Long>();
        nLoci = 0;
        totalLocusDepth = 0;
        totalDepthOfCoverage = 0;
    }

    public void addSample(String sample) {
        if ( granularHistogramBySample.containsKey(sample) ) {
            return;
        }

        int[] binCounts = new int[this.binLeftEndpoints.length+1];
        for ( int b = 0; b < binCounts.length; b ++ ) {
            binCounts[b] = 0;
        }

        granularHistogramBySample.put(sample,binCounts);
        totalCoverages.put(sample,0l);
    }

    public void initializeLocusCounts() {
        locusCoverageCounts = new int[granularHistogramBySample.size()][binLeftEndpoints.length+1];
        locusHistogram = new int[binLeftEndpoints.length+1];
        for ( int b = 0; b < binLeftEndpoints.length+1; b ++ ) {
            for ( int a = 0; a < granularHistogramBySample.size(); a ++ ) {
                locusCoverageCounts[a][b] = 0;
            }
            locusHistogram[b] = 0;
        }

        tabulateLocusCounts = true;
    }

    public void initializeDeletions() {
        includeDeletions = true;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // UPDATE METHODS
    ////////////////////////////////////////////////////////////////////////////////////

    public void updateDepths(Map<String,Integer> depthBySample) {
        int b;
        for ( String sample : granularHistogramBySample.keySet() ) {
            if ( depthBySample.containsKey(sample) ) {
                b = updateSample(sample,depthBySample.get(sample));
                totalLocusDepth += depthBySample.get(sample);
            } else {
                b = updateSample(sample,0);
            }

            if ( tabulateLocusCounts ) {
                for ( int i = 0; i <= b; i ++ ) {
                    locusHistogram[i]++;
                }
            }
        }
        updateLocusCounts(locusHistogram);

        nLoci++;
        totalDepthOfCoverage += totalLocusDepth;
        totalLocusDepth = 0;
    }

    public void update(Map<String,int[]> countsBySample) {
        // todo -- do we want to do anything special regarding base count or deletion statistics?
        HashMap<String,Integer> depthBySample = new HashMap<String,Integer>();
        // todo -- needs fixing with advent of new baseutils functionality using ENUMS and handling N,D
        for ( String s : countsBySample.keySet() ) {
            int total = 0;
            int[] counts = countsBySample.get(s);
            for ( char base : BaseUtils.EXTENDED_BASES ) {
                if ( includeDeletions || ! ( base == 'D') ) { // note basesAreEqual assigns TRUE to (N,D) as both have simple index -1
                    total += counts[BaseUtils.extendedBaseToBaseIndex(base)];
                }
            }
            depthBySample.put(s,total);
        }
        
        this.updateDepths(depthBySample);
    }

    private int updateSample(String sample, int depth) {
        totalCoverages.put(sample,totalCoverages.get(sample)+depth);

        int[] granularBins = granularHistogramBySample.get(sample);
        for ( int b = 0; b < binLeftEndpoints.length; b ++ ) {
            if ( depth < binLeftEndpoints[b] ) {
                granularBins[b]++;
                return b;
            }
        }

        granularBins[binLeftEndpoints.length]++; // greater than all left-endpoints
        return binLeftEndpoints.length;
    }

    public void merge(DepthOfCoverageStats newStats) {
        this.mergeSamples(newStats);
        if ( this.tabulateLocusCounts && newStats.tabulateLocusCounts ) {
            this.mergeLocusCounts(newStats.getLocusCounts());
        }
        nLoci += newStats.getTotalLoci();
        totalDepthOfCoverage += newStats.getTotalCoverage();
    }

    private void mergeSamples(DepthOfCoverageStats otherStats) {
        Map<String,int[]> otherHistogram = otherStats.getHistograms();
        Map<String,Double> otherMeans = otherStats.getMeans();
        for ( String s : this.getAllSamples() ) {
            int[] internalCounts = granularHistogramBySample.get(s);
            int[] externalCounts = otherHistogram.get(s);
            for ( int b = 0; b < internalCounts.length; b++ ) {
                internalCounts[b] += externalCounts[b];
            }

            this.totalCoverages.put(s, this.totalCoverages.get(s) + otherStats.totalCoverages.get(s));
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
        HashMap<String,Double> means = new HashMap<String,Double>();
        for ( String s : getAllSamples() ) {
            means.put(s,( (double)totalCoverages.get(s))/( (double) nLoci ));
        }

        return means;
    }

    public Map<String,Long> getTotals() {
        return totalCoverages;
    }

    public long getTotalLoci() {
        return nLoci;
    }

    public Set<String> getAllSamples() {
        return granularHistogramBySample.keySet();
    }

    public double getTotalMeanCoverage() {
        return ( (double) totalDepthOfCoverage )/ ( (double) nLoci );
    }

    public long getTotalCoverage() {
        return totalDepthOfCoverage;
    }

}