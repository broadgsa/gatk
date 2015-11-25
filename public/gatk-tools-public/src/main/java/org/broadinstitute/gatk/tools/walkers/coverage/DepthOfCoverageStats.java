/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

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

    private Map<String,long[]> granularHistogramBySample; // holds the counts per each bin
    private Map<String,Long> totalCoverages; // holds total coverage per sample
    private int[] binLeftEndpoints; // describes the left endpoint for each bin
    private long[][] locusCoverageCounts; // holds counts of number of bases with >=X samples at >=Y coverage
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
            throw new UserException.BadInput("the start must be at least 1 and the number of bins may not exceed stop - start");
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
        granularHistogramBySample = new HashMap<String,long[]>();
        totalCoverages = new HashMap<String,Long>();
        nLoci = 0;
        totalLocusDepth = 0;
        totalDepthOfCoverage = 0;
    }

    public DepthOfCoverageStats(DepthOfCoverageStats cloneMe) {
        this.binLeftEndpoints = cloneMe.binLeftEndpoints;
        granularHistogramBySample = new TreeMap<String,long[]>();
        totalCoverages = new TreeMap<String,Long>();
        for ( String s : cloneMe.getAllSamples() ) {
            granularHistogramBySample.put(s,new long[cloneMe.getHistograms().get(s).length]);
            for ( int i = 0; i < granularHistogramBySample.get(s).length; i++ ) {
                granularHistogramBySample.get(s)[i] = cloneMe.getHistograms().get(s)[i];
            }
            totalCoverages.put(s,cloneMe.totalCoverages.get(s));
        }

        this.includeDeletions = cloneMe.includeDeletions;
        if ( cloneMe.tabulateLocusCounts ) {
            this.locusCoverageCounts = new long[cloneMe.locusCoverageCounts.length][cloneMe.locusCoverageCounts[0].length];
        }
        //this.granularHistogramBySample = cloneMe.granularHistogramBySample;
        //this.totalCoverages = cloneMe.totalCoverages;
        this.nLoci = cloneMe.nLoci;
        this.totalDepthOfCoverage = cloneMe.totalDepthOfCoverage;
        this.tabulateLocusCounts = cloneMe.tabulateLocusCounts;
    }

    public void addSample(String sample) {
        if ( granularHistogramBySample.containsKey(sample) ) {
            return;
        }

        long[] binCounts = new long[this.binLeftEndpoints.length+1];
        for ( int b = 0; b < binCounts.length; b ++ ) {
            binCounts[b] = 0;
        }

        granularHistogramBySample.put(sample,binCounts);
        totalCoverages.put(sample,0l);
    }

    public void initializeLocusCounts() {
        locusCoverageCounts = new long[granularHistogramBySample.size()][binLeftEndpoints.length+1];
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
        if ( countsBySample == null ) {
            this.updateDepths(new HashMap<String,Integer>(1));
            return;
        }
        // todo -- do we want to do anything special regarding base count or deletion statistics?
        HashMap<String,Integer> depthBySample = new HashMap<String,Integer>();
        // todo -- needs fixing with advent of new baseutils functionality using ENUMS and handling N,D
        for ( String s : countsBySample.keySet() ) {
            int total = 0;
            int[] counts = countsBySample.get(s);
            for ( byte base : BaseUtils.EXTENDED_BASES ) {
                if ( includeDeletions || ! ( base == BaseUtils.Base.D.base) ) { // note basesAreEqual assigns TRUE to (N,D) as both have simple index -1
                    total += counts[BaseUtils.extendedBaseToBaseIndex(base)];
                }
            }
            depthBySample.put(s,total);
        }
        
        this.updateDepths(depthBySample);
    }

    private int updateSample(String sample, int depth) {
        totalCoverages.put(sample,totalCoverages.get(sample)+depth);

        long[] granularBins = granularHistogramBySample.get(sample);
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
        Map<String,long[]> otherHistogram = otherStats.getHistograms();
        Map<String,Double> otherMeans = otherStats.getMeans();
        for ( String s : this.getAllSamples() ) {
            long[] internalCounts = granularHistogramBySample.get(s);
            long[] externalCounts = otherHistogram.get(s);
            for ( int b = 0; b < internalCounts.length; b++ ) {
                internalCounts[b] += externalCounts[b];
            }

            this.totalCoverages.put(s, this.totalCoverages.get(s) + otherStats.totalCoverages.get(s));
        }
    }

    private void mergeLocusCounts( long[][] otherCounts ) {
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

    public Map<String,long[]> getHistograms() {
        return granularHistogramBySample;
    }

    public long[][] getLocusCounts() {
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

    public double[] getCoverageProportions(String sample) {
        long[] hist = granularHistogramBySample.get(sample);
        double[] distribution = new double[hist.length];
        long count = 0;
        for ( int i = hist.length-1; i >= 0; i -- ) {
            count += hist[i];
            distribution[i] = ( (double) count) / nLoci;
        }

        return distribution;
    }

    public int value2bin(int value) {
        for ( int index = 0; index < binLeftEndpoints.length; index++ ) {
            if ( binLeftEndpoints[index] > value ) {
                return index;
            }
        }

        return binLeftEndpoints.length-1;
    }

}