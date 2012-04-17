/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.*;

/**
 * Util class for performing the pair HMM for local alignment. Figure 4.3 in Durbin 1998 book.
 * User: rpoplin
 * Date: 3/1/12
 */

public class PairHMM {
    private static final int MAX_CACHED_QUAL = (int)Byte.MAX_VALUE;
    private static final byte DEFAULT_GOP = (byte) 45;
    private static final byte DEFAULT_GCP = (byte) 10;
    private static final double BANDING_TOLERANCE = 22.0;
    private static final int BANDING_CLUSTER_WINDOW = 12;
    private final boolean doBanded;

    public PairHMM() {
        doBanded = false;
    }

    public PairHMM( final boolean doBanded ) {
        this.doBanded = doBanded;
    }

    
    public static void initializeArrays(final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray,
                                 final int X_METRIC_LENGTH) {

        for( int iii=0; iii < X_METRIC_LENGTH; iii++ ) {
            Arrays.fill(matchMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[iii], Double.NEGATIVE_INFINITY);
        }

        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);

    }

    @Requires({"readBases.length == readQuals.length","readBases.length == insertionGOP.length","readBases.length == deletionGOP.length","readBases.length == overallGCP.length"})
    @Ensures({"!Double.isInfinite(result)", "!Double.isNaN(result)"}) // Result should be a proper log10 probability
    public double computeReadLikelihoodGivenHaplotype( final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals,
                                                       final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP ) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions
        final int X_METRIC_LENGTH = readBases.length + 1;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 1;

        // initial arrays to hold the probabilities of being in the match, insertion and deletion cases
        final double[][] matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        final double[][] YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        initializeArrays(matchMetricArray, XMetricArray, YMetricArray, X_METRIC_LENGTH);

        return computeReadLikelihoodGivenHaplotype(haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 0, matchMetricArray, XMetricArray, YMetricArray);
    }

    @Requires({"readBases.length == readQuals.length","readBases.length == insertionGOP.length","readBases.length == deletionGOP.length","readBases.length == overallGCP.length"})
    @Ensures({"!Double.isInfinite(result)", "!Double.isNaN(result)"}) // Result should be a proper log10 probability
    public double computeReadLikelihoodGivenHaplotype( final byte[] haplotypeBases, final byte[] readBases, final byte[] readQuals,
                                                       final byte[] insertionGOP, final byte[] deletionGOP, final byte[] overallGCP, final int hapStartIndex,
                                                       final double[][] matchMetricArray, final double[][] XMetricArray, final double[][] YMetricArray ) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions
        final int X_METRIC_LENGTH = readBases.length + 1;
        final int Y_METRIC_LENGTH = haplotypeBases.length + 1;

        if( doBanded ) {
            final ArrayList<Integer> workQueue = new ArrayList<Integer>(); // holds a queue of starting work location (indices along the diagonal). Will be sorted each step
            final ArrayList<Integer> workToBeAdded = new ArrayList<Integer>();
            final ArrayList<Double> calculatedValues = new ArrayList<Double>();
            final int numDiags = X_METRIC_LENGTH + Y_METRIC_LENGTH - 1;
            workQueue.add( 1 ); // Always start a new thread at the baseline because of partially repeating sequences that match better in the latter half of the haplotype

            for(int diag = 3; diag < numDiags; diag++) { // diag = 3 is the (1,2) element of the metric arrays. (1,1) is the initial condition and is purposefully skipped over
                //Collections.sort(workQueue); // no need to sort because elements are guaranteed to be in ascending order
                int el = 1;
                for( int work : workQueue ) {
                    // choose the appropriate diagonal baseline location
                    int iii = 0;
                    int jjj = diag;
                    if( diag > Y_METRIC_LENGTH ) {
                        iii = diag - Y_METRIC_LENGTH;
                        jjj = Y_METRIC_LENGTH;
                    }
                    // move to the starting work location along the diagonal
                    iii += work;
                    jjj -= work;
                    while( iii >= X_METRIC_LENGTH || jjj <= 0 ) {
                        iii--;
                        jjj++;
                        work--;
                    }
                    if( !detectClusteredStartLocations(workToBeAdded, work ) ) {
                        workToBeAdded.add(work); // keep this thread going once it has started
                    }

                    if( work >= el - 3 ) {
                        // step along the diagonal in the forward direction, updating the match matrices and looking for a drop off from the maximum observed value
                        double maxElement = Double.NEGATIVE_INFINITY;
                        for( el = work; el < numDiags + 1; el++ ) {
                            updateCell(iii, jjj, haplotypeBases, readBases, readQuals,
                                insertionGOP, deletionGOP, overallGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                            final double bestMetric = MathUtils.max(matchMetricArray[iii][jjj], XMetricArray[iii][jjj], YMetricArray[iii][jjj]);
                            calculatedValues.add(bestMetric);
                            if( bestMetric > maxElement ) {
                                maxElement = bestMetric;
                            } else if( maxElement - bestMetric > BANDING_TOLERANCE ) {
                                break;
                            }
                            if( ++iii >= X_METRIC_LENGTH ) { // don't walk off the edge of the matrix
                                break;
                            }
                            if( --jjj <= 0 ) { // don't walk off the edge of the matrix
                                break;
                            }
                        }

                        // find a local maximum to start a new band in the work queue
                        double localMaxElement = Double.NEGATIVE_INFINITY;
                        int localMaxElementIndex = 0;
                        for(int kkk = calculatedValues.size()-1; kkk >= 1; kkk--) {
                            final double bestMetric = calculatedValues.get(kkk);
                            if( bestMetric > localMaxElement ) {
                                localMaxElement = bestMetric;
                                localMaxElementIndex = kkk;
                            } else if( localMaxElement - bestMetric > BANDING_TOLERANCE * 0.5 ) { // find a local maximum
                                if( !detectClusteredStartLocations(workToBeAdded, work + localMaxElementIndex ) ) {
                                    workToBeAdded.add( work + localMaxElementIndex );
                                }
                                break;
                            }
                        }
                        calculatedValues.clear();

                        // reset iii and jjj to the appropriate diagonal baseline location
                        iii = 0;
                        jjj = diag;
                        if( diag > Y_METRIC_LENGTH ) {
                            iii = diag - Y_METRIC_LENGTH;
                            jjj = Y_METRIC_LENGTH;
                        }
                        // move to the starting work location along the diagonal
                        iii += work-1;
                        jjj -= work-1;

                        // step along the diagonal in the reverse direction, updating the match matrices and looking for a drop off from the maximum observed value
                        for( int traceBack = work - 1; traceBack > 0 && iii > 0 && jjj < Y_METRIC_LENGTH; traceBack--,iii--,jjj++ ) {
                            updateCell(iii, jjj, haplotypeBases, readBases, readQuals,
                                    insertionGOP, deletionGOP, overallGCP,  matchMetricArray,  XMetricArray, YMetricArray);
                            final double bestMetric = MathUtils.max(matchMetricArray[iii][jjj], XMetricArray[iii][jjj], YMetricArray[iii][jjj]);
                            if( bestMetric > maxElement ) {
                                maxElement = bestMetric;
                            } else if( maxElement - bestMetric > BANDING_TOLERANCE ) {
                                break;
                            }
                        }
                    }
                }
                workQueue.clear();
                workQueue.addAll(workToBeAdded);
                workToBeAdded.clear();
            }
        } else {
            // simple rectangular version of update loop, slow
            for( int iii = 1; iii < X_METRIC_LENGTH; iii++ ) {
                for( int jjj = hapStartIndex + 1; jjj < Y_METRIC_LENGTH; jjj++ ) {
                    if( (iii == 1 && jjj == 1) ) { continue; }
                    updateCell(iii, jjj, haplotypeBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP,
                        matchMetricArray, XMetricArray, YMetricArray);
                }
            }
        }

        // final probability is the log10 sum of the last element in all three state arrays
        final int endI = X_METRIC_LENGTH - 1;
        final int endJ = Y_METRIC_LENGTH - 1;
        return MathUtils.approximateLog10SumLog10(new double[]{matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]});
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
            final byte qual = ( readQuals[im1-1] < QualityUtils.MIN_USABLE_Q_SCORE ? QualityUtils.MIN_USABLE_Q_SCORE : (readQuals[im1-1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : readQuals[im1-1]) );
            pBaseReadLog10 = ( x == y || x == (byte) 'N' || y == (byte) 'N' ? QualityUtils.qualToProbLog10(qual) : QualityUtils.qualToErrorProbLog10(qual) );
        }
        final int qualIndexGOP = ( im1 == 0 ? DEFAULT_GOP + DEFAULT_GOP : ( insertionGOP[im1-1] + deletionGOP[im1-1] > MAX_CACHED_QUAL ? MAX_CACHED_QUAL : insertionGOP[im1-1] + deletionGOP[im1-1]) );
        final double d0 = QualityUtils.qualToProbLog10((byte)qualIndexGOP);
        final double e0 = ( im1 == 0 ? QualityUtils.qualToProbLog10(DEFAULT_GCP) : QualityUtils.qualToProbLog10(overallGCP[im1-1]) );
        matchMetricArray[indI][indJ] = pBaseReadLog10 + MathUtils.approximateLog10SumLog10(
                new double[]{matchMetricArray[indI-1][indJ-1] + d0, XMetricArray[indI-1][indJ-1] + e0, YMetricArray[indI-1][indJ-1] + e0});

        // update the X (insertion) array
        final double d1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GOP) : QualityUtils.qualToErrorProbLog10(insertionGOP[im1-1]) );
        final double e1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GCP) : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseReadLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        XMetricArray[indI][indJ] = qBaseReadLog10 + MathUtils.approximateLog10SumLog10(matchMetricArray[indI-1][indJ] + d1, XMetricArray[indI-1][indJ] + e1);

        // update the Y (deletion) array, with penalty of zero on the left and right flanks to allow for a local alignment within the haplotype
        final double d2 = ( im1 == 0 || im1 == readBases.length - 1 ? 0.0 : QualityUtils.qualToErrorProbLog10(deletionGOP[im1-1]) );
        final double e2 = ( im1 == 0 || im1 == readBases.length - 1 ? 0.0 : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseRefLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        YMetricArray[indI][indJ] = qBaseRefLog10 + MathUtils.approximateLog10SumLog10(matchMetricArray[indI][indJ-1] + d2, YMetricArray[indI][indJ-1] + e2);
    }

    // private function used by the banded approach to ensure the proposed bands are sufficiently distinct from each other
    private boolean detectClusteredStartLocations( final ArrayList<Integer> list, int loc ) {
        for(int x : list) {
            if( Math.abs(x-loc) <= BANDING_CLUSTER_WINDOW ) {
                return true;
            }
        }
        return false;
    }
}