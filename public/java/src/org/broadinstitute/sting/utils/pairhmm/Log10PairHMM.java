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

import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.Arrays;

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 *
 * User: rpoplin
 * Date: 3/1/12
 */
public class Log10PairHMM extends PairHMM {
    /**
     * Should we use exact log10 calculation (true), or an approximation (false)?
     */
    private final boolean doExactLog10;

    /**
     * Create an uninitialized PairHMM
     *
     * @param doExactLog10 should the log10 calculations be exact (slow) or approximate (faster)
     */
    public Log10PairHMM(final boolean doExactLog10) {
        this.doExactLog10 = doExactLog10;
    }

    /**
     * Is this HMM using exact log10 calculations?
     * @return true if exact, false if approximate
     */
    public boolean isDoingExactLog10Calculations() {
        return doExactLog10;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize( final int readMaxLength, final int haplotypeMaxLength) {
        super.initialize(readMaxLength, haplotypeMaxLength);

        for( int iii=0; iii < X_METRIC_MAX_LENGTH; iii++ ) {
            Arrays.fill(matchMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[iii], Double.NEGATIVE_INFINITY);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                               final byte[] readBases,
                                                               final byte[] readQuals,
                                                               final byte[] insertionGOP,
                                                               final byte[] deletionGOP,
                                                               final byte[] overallGCP,
                                                               final int hapStartIndex,
                                                               final boolean recacheReadValues ) {
        // the initial condition -- must be in subComputeReadLikelihoodGivenHaplotypeLog10 because it needs that actual
        // read and haplotypes, not the maximum
        matchMetricArray[1][1] = getNPotentialXStartsLikelihoodPenaltyLog10(haplotypeBases.length, readBases.length);

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = readBases.length + 2;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 2;

        // ensure that all the qual scores have valid values
        for( int iii = 0; iii < readQuals.length; iii++ ) {
            readQuals[iii] = ( readQuals[iii] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[iii] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[iii]) );
        }

        // simple rectangular version of update loop, slow
        for( int iii = 1; iii < X_METRIC_LENGTH; iii++ ) {
            for( int jjj = hapStartIndex + 1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                if( (iii == 1 && jjj == 1) ) { continue; }
                updateCell(iii, jjj, haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP,
                        matchMetricArray, XMetricArray, YMetricArray);
            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        return myLog10SumLog10(new double[]{matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]});
    }

    /**
     * Compute the log10SumLog10 of the values
     *
     * NOTE NOTE NOTE
     *
     * Log10PairHMM depends critically on this function tolerating values that are all -Infinity
     * and the sum returning -Infinity.  Note good.  Needs to be fixed.
     *
     * NOTE NOTE NOTE
     *
     * @param values an array of log10 probabilities that need to be summed
     * @return the log10 of the sum of the probabilities
     */
    @Requires("values != null")
    private double myLog10SumLog10(final double[] values) {
        return doExactLog10 ? MathUtils.log10sumLog10(values) : MathUtils.approximateLog10SumLog10(values);
    }

    private void updateCell( final int indI, final int indJ, final byte[] haplotypeBases, final byte[] readBases,
                             final byte[] readQuals, final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP,
                             final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        // the read and haplotype indices are offset by one because the state arrays have an extra column to hold the initial conditions
        final int im1 = indI - 1;
        final int jm1 = indJ - 1;

        // update the match array
        double pBaseReadLog10 = 0.0; // Math.log10(1.0);
        if( im1 > 0 && jm1 > 0 ) { // the emission probability is applied when leaving the state
            final byte x = readBases[im1-1];
            final byte y = haplotypeBases[jm1-1];
            final byte qual = readQuals[im1-1];
            pBaseReadLog10 = ( x == y || x == (byte) 'N' || y == (byte) 'N' ? QualityUtils.qualToProbLog10(qual) : QualityUtils.qualToErrorProbLog10(qual) );
        }
        final int qualIndexGOP = ( im1 == 0 ? DEFAULT_GOP + DEFAULT_GOP : ( insertionGOP[im1-1] + deletionGOP[im1-1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : insertionGOP[im1-1] + deletionGOP[im1-1]) );
        final double d0 = QualityUtils.qualToProbLog10((byte)qualIndexGOP);
        final double e0 = ( im1 == 0 ? QualityUtils.qualToProbLog10(DEFAULT_GCP) : QualityUtils.qualToProbLog10(overallGCP[im1-1]) );
        matchMetricArray[indI][indJ] = pBaseReadLog10 + myLog10SumLog10(new double[]{matchMetricArray[indI - 1][indJ - 1] + d0, XMetricArray[indI - 1][indJ - 1] + e0, YMetricArray[indI - 1][indJ - 1] + e0});

        // update the X (insertion) array
        final double d1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GOP) : QualityUtils.qualToErrorProbLog10(insertionGOP[im1-1]) );
        final double e1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GCP) : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseReadLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        XMetricArray[indI][indJ] = qBaseReadLog10 + myLog10SumLog10(new double[]{matchMetricArray[indI - 1][indJ] + d1, XMetricArray[indI - 1][indJ] + e1});

        // update the Y (deletion) array, with penalty of zero on the left and right flanks to allow for a local alignment within the haplotype
        final double d2 = ( im1 == 0 || im1 == readBases.length ? 0.0 : QualityUtils.qualToErrorProbLog10(deletionGOP[im1-1]) );
        final double e2 = ( im1 == 0 || im1 == readBases.length ? 0.0 : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseRefLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        YMetricArray[indI][indJ] = qBaseRefLog10 + myLog10SumLog10(new double[]{matchMetricArray[indI][indJ - 1] + d2, YMetricArray[indI][indJ - 1] + e2});
    }
}
