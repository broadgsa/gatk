/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

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
        granularHistogramBySample = new HashMap<String,long[]>();
        totalCoverages = new HashMap<String,Long>();
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