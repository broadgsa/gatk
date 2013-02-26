/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.pairhmm;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 *
 * User: rpoplin
 * Date: 10/16/12
 */
public abstract class PairHMM {
    protected final static Logger logger = Logger.getLogger(PairHMM.class);

    protected static final Byte MAX_CACHED_QUAL = Byte.MAX_VALUE;
    protected static final byte DEFAULT_GOP = (byte) 45;
    protected static final byte DEFAULT_GCP = (byte) 10;

    public enum HMM_IMPLEMENTATION {
        /* Very slow implementation which uses very accurate log10 sum functions. Only meant to be used as a reference test implementation */
        EXACT,
        /* PairHMM as implemented for the UnifiedGenotyper. Uses log10 sum functions accurate to only 1E-4 */
        ORIGINAL,
        /* Optimized version of the PairHMM which caches per-read computations and operations in real space to avoid costly sums of log10'ed likelihoods */
        LOGLESS_CACHING
    }

    protected double[][] matchMetricArray = null;
    protected double[][] XMetricArray = null;
    protected double[][] YMetricArray = null;
    protected int maxHaplotypeLength, maxReadLength;
    protected int X_METRIC_MAX_LENGTH, Y_METRIC_MAX_LENGTH;
    private boolean initialized = false;

    /**
     * Initialize this PairHMM, making it suitable to run against a read and haplotype with given lengths
     * @param readMaxLength the max length of reads we want to use with this PairHMM
     * @param haplotypeMaxLength the max length of haplotypes we want to use with this PairHMM
     */
    public void initialize( final int readMaxLength, final int haplotypeMaxLength ) {
        if ( readMaxLength <= 0 ) throw new IllegalArgumentException("READ_MAX_LENGTH must be > 0 but got " + readMaxLength);
        if ( haplotypeMaxLength <= 0 ) throw new IllegalArgumentException("HAPLOTYPE_MAX_LENGTH must be > 0 but got " + haplotypeMaxLength);

        maxHaplotypeLength = haplotypeMaxLength;
        maxReadLength = readMaxLength;

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        X_METRIC_MAX_LENGTH = readMaxLength + 2;
        Y_METRIC_MAX_LENGTH = haplotypeMaxLength + 2;

        matchMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        XMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        YMetricArray = new double[X_METRIC_MAX_LENGTH][Y_METRIC_MAX_LENGTH];
        initialized = true;
    }

    /**
     * Compute the total probability of read arising from haplotypeBases given base substitution, insertion, and deletion
     * probabilities.
     *
     * Note on using hapStartIndex.  This allows you to compute the exact true likelihood of a full haplotypes
     * given a read, assuming that the previous calculation read over a full haplotype, recaching the read values,
     * starting only at the place where the new haplotype bases and the previous haplotype bases different.  This
     * index is 0-based, and can be computed with findFirstPositionWhereHaplotypesDiffer given the two haplotypes.
     * Note that this assumes that the read and all associated quals values are the same.
     *
     * @param haplotypeBases the full sequence (in standard SAM encoding) of the haplotype, must be >= than read bases in length
     * @param readBases the bases (in standard encoding) of the read, must be <= haplotype bases in length
     * @param readQuals the phred-scaled per base substitition quality scores of read.  Must be the same length as readBases
     * @param insertionGOP the phred-scaled per base insertion quality scores of read.  Must be the same length as readBases
     * @param deletionGOP the phred-scaled per base deletion quality scores of read.  Must be the same length as readBases
     * @param overallGCP the phred-scaled gap continuation penalties scores of read.  Must be the same length as readBases
     * @param hapStartIndex start the hmm calculation at this offset in haplotype bases.  Used in the caching calculation
     *                      where multiple haplotypes are used, and they only diff starting at hapStartIndex
     * @param recacheReadValues if false, we don't recalculate any cached results, assuming that readBases and its associated
     *                          parameters are the same, and only the haplotype bases are changing underneath us
     * @return the log10 probability of read coming from the haplotype under the provided error model
     */
    public final double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                  final byte[] readBases,
                                                                  final byte[] readQuals,
                                                                  final byte[] insertionGOP,
                                                                  final byte[] deletionGOP,
                                                                  final byte[] overallGCP,
                                                                  final int hapStartIndex,
                                                                  final boolean recacheReadValues ) {
        if ( ! initialized ) throw new IllegalStateException("Must call initialize before calling computeReadLikelihoodGivenHaplotypeLog10");
        if ( haplotypeBases == null ) throw new IllegalArgumentException("haplotypeBases cannot be null");
        if ( haplotypeBases.length > maxHaplotypeLength ) throw new IllegalArgumentException("Haplotype bases is too long, got " + haplotypeBases.length + " but max is " + maxHaplotypeLength);
        if ( readBases == null ) throw new IllegalArgumentException("readBases cannot be null");
        if ( readBases.length > maxReadLength ) throw new IllegalArgumentException("readBases is too long, got " + readBases.length + " but max is " + maxReadLength);
        if ( readQuals.length != readBases.length ) throw new IllegalArgumentException("Read bases and read quals aren't the same size: " + readBases.length + " vs " + readQuals.length);
        if ( insertionGOP.length != readBases.length ) throw new IllegalArgumentException("Read bases and read insertion quals aren't the same size: " + readBases.length + " vs " + insertionGOP.length);
        if ( deletionGOP.length != readBases.length ) throw new IllegalArgumentException("Read bases and read deletion quals aren't the same size: " + readBases.length + " vs " + deletionGOP.length);
        if ( overallGCP.length != readBases.length ) throw new IllegalArgumentException("Read bases and overall GCP aren't the same size: " + readBases.length + " vs " + overallGCP.length);
        if ( hapStartIndex < 0 || hapStartIndex > haplotypeBases.length ) throw new IllegalArgumentException("hapStartIndex is bad, must be between 0 and haplotype length " + haplotypeBases.length + " but got " + hapStartIndex);

        final double result = subComputeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, hapStartIndex, recacheReadValues);

        if ( MathUtils.goodLog10Probability(result) )
            return result;
        else
            throw new IllegalStateException("Bad likelihoods detected: " + result);
//            return result;
    }

    /**
     * To be overloaded by subclasses to actually do calculation for #computeReadLikelihoodGivenHaplotypeLog10
     */
    @Requires({"readBases.length == readQuals.length", "readBases.length == insertionGOP.length", "readBases.length == deletionGOP.length",
               "readBases.length == overallGCP.length", "matchMetricArray!=null", "XMetricArray!=null", "YMetricArray!=null"})
    protected abstract double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                        final byte[] readBases,
                                                                        final byte[] readQuals,
                                                                        final byte[] insertionGOP,
                                                                        final byte[] deletionGOP,
                                                                        final byte[] overallGCP,
                                                                        final int hapStartIndex,
                                                                        final boolean recacheReadValues );

    /**
     * How many potential starting locations are a read with readSize bases against a haplotype with haplotypeSize bases?
     *
     * for example, a 3 bp read against a 5 bp haplotype could potentially start at 1, 2, 3 = 5 - 3 + 1 = 3
     * the max value is necessary in the case where the read is longer than the haplotype, in which case
     * there's a single unique start site by assumption
     *
     * @param haplotypeSize the number of bases in the haplotype we are testing
     * @param readSize the number of bases in the read we are testing
     * @return a positive integer >= 1
     */
    @Ensures("result >= 1")
    protected int getNPotentialXStarts(final int haplotypeSize, final int readSize) {
        return Math.max(haplotypeSize - readSize + 1, 1);
    }

    /**
     * The the log10 probability penalty for the number of potential start sites of the read aginst the haplotype
     *
     * @param haplotypeSize the number of bases in the haplotype we are testing
     * @param readSize the number of bases in the read we are testing
     * @return a log10 probability
     */
    @Ensures("MathUtils.goodLog10Probability(result)")
    protected double getNPotentialXStartsLikelihoodPenaltyLog10(final int haplotypeSize, final int readSize) {
        return - Math.log10(getNPotentialXStarts(haplotypeSize, readSize));
    }

    /**
     * Print out the core hmm matrices for debugging
     */
    protected void dumpMatrices() {
        dumpMatrix("matchMetricArray", matchMetricArray);
        dumpMatrix("XMetricArray", XMetricArray);
        dumpMatrix("YMetricArray", YMetricArray);
    }

    /**
     * Print out in a human readable form the matrix for debugging
     * @param name the name of this matrix
     * @param matrix the matrix of values
     */
    @Requires({"name != null", "matrix != null"})
    private void dumpMatrix(final String name, final double[][] matrix) {
        System.out.printf("%s%n", name);
        for ( int i = 0; i < matrix.length; i++) {
            System.out.printf("\t%s[%d]", name, i);
            for ( int j = 0; j < matrix[i].length; j++ ) {
                if ( Double.isInfinite(matrix[i][j]) )
                    System.out.printf(" %15s", String.format("%f", matrix[i][j]));
                else
                    System.out.printf(" % 15.5e", matrix[i][j]);
            }
            System.out.println();
        }
    }

    /**
     * Compute the first position at which two haplotypes differ
     *
     * If the haplotypes are exact copies of each other, returns the min length of the two haplotypes.
     *
     * @param haplotype1 the first haplotype1
     * @param haplotype2 the second haplotype1
     * @return the index of the first position in haplotype1 and haplotype2 where the byte isn't the same
     */
    public static int findFirstPositionWhereHaplotypesDiffer(final byte[] haplotype1, final byte[] haplotype2) {
        if ( haplotype1 == null || haplotype1.length == 0 ) throw new IllegalArgumentException("Haplotype1 is bad " + haplotype1);
        if ( haplotype2 == null || haplotype2.length == 0 ) throw new IllegalArgumentException("Haplotype2 is bad " + haplotype2);

        for( int iii = 0; iii < haplotype1.length && iii < haplotype2.length; iii++ ) {
            if( haplotype1[iii] != haplotype2[iii] ) {
                return iii;
            }
        }

        return Math.min(haplotype1.length, haplotype2.length);
    }
}
