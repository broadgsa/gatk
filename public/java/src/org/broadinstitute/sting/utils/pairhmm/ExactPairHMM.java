package org.broadinstitute.sting.utils.pairhmm;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 10/16/12
 */

public class ExactPairHMM extends PairHMM {

    @Override
    public void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH ) {

        // M, X, and Y arrays are of size read and haplotype + 1 because of an extra column for initial conditions and + 1 to consider the final base in a non-global alignment
        final int X_METRIC_LENGTH = READ_MAX_LENGTH + 2;
        final int Y_METRIC_LENGTH = HAPLOTYPE_MAX_LENGTH + 2;

        matchMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        XMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];
        YMetricArray = new double[X_METRIC_LENGTH][Y_METRIC_LENGTH];

        for( int iii=0; iii < X_METRIC_LENGTH; iii++ ) {
            Arrays.fill(matchMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(XMetricArray[iii], Double.NEGATIVE_INFINITY);
            Arrays.fill(YMetricArray[iii], Double.NEGATIVE_INFINITY);
        }

        // the initial condition
        matchMetricArray[1][1] = 0.0; // Math.log10(1.0);
    }

    @Override
    public double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                            final byte[] readBases,
                                                            final byte[] readQuals,
                                                            final byte[] insertionGOP,
                                                            final byte[] deletionGOP,
                                                            final byte[] overallGCP,
                                                            final int hapStartIndex,
                                                            final boolean recacheReadValues ) {

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
        return MathUtils.log10sumLog10(new double[]{matchMetricArray[endI][endJ], XMetricArray[endI][endJ], YMetricArray[endI][endJ]});
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
        matchMetricArray[indI][indJ] = pBaseReadLog10 + MathUtils.log10sumLog10(new double[]{matchMetricArray[indI-1][indJ-1] + d0, XMetricArray[indI-1][indJ-1] + e0, YMetricArray[indI-1][indJ-1] + e0});

        // update the X (insertion) array
        final double d1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GOP) : QualityUtils.qualToErrorProbLog10(insertionGOP[im1-1]) );
        final double e1 = ( im1 == 0 ? QualityUtils.qualToErrorProbLog10(DEFAULT_GCP) : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseReadLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        XMetricArray[indI][indJ] = qBaseReadLog10 + MathUtils.log10sumLog10(new double[]{matchMetricArray[indI-1][indJ] + d1, XMetricArray[indI-1][indJ] + e1});

        // update the Y (deletion) array, with penalty of zero on the left and right flanks to allow for a local alignment within the haplotype
        final double d2 = ( im1 == 0 || im1 == readBases.length ? 0.0 : QualityUtils.qualToErrorProbLog10(deletionGOP[im1-1]) );
        final double e2 = ( im1 == 0 || im1 == readBases.length ? 0.0 : QualityUtils.qualToErrorProbLog10(overallGCP[im1-1]) );
        final double qBaseRefLog10 = 0.0; // Math.log10(1.0) -- we don't have an estimate for this emission probability so assume q=1.0
        YMetricArray[indI][indJ] = qBaseRefLog10 + MathUtils.log10sumLog10(new double[]{matchMetricArray[indI][indJ-1] + d2, YMetricArray[indI][indJ-1] + e2});
    }
}
