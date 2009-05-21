package org.broadinstitute.sting.secondarybase;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * BasecallingBaseModel is a class that represents the statistical
 * model for all bases at a given cycle.  It allows for easy, one
 * pass training via the addTrainingPoint() method.  Once the model
 * is trained, computeLikelihoods will return the probability matrix
 * over previous cycle's base hypotheses and current cycle base
 * hypotheses (contextual prior is included in these likelihoods).
 *
 * @author Kiran Garimella
 */
public class BasecallingBaseModel {
    private double[][] counts;
    private DoubleMatrix1D[][] sums;
    private DoubleMatrix2D[][] unscaledCovarianceSums;

    private DoubleMatrix1D[][] means;
    private DoubleMatrix2D[][] inverseCovariances;
    private double[][] norms;

    private cern.jet.math.Functions F = cern.jet.math.Functions.functions;
    private Algebra alg;

    private boolean correctForContext = false;
    private int numTheories = 1;

    private boolean readyToCall = false;

    /**
     * Constructor for BasecallingBaseModel.
     *
     * @param correctForContext  should we attempt to correct for contextual sequence effects?
     */
    public BasecallingBaseModel(boolean correctForContext) {
        this.correctForContext = correctForContext;
        this.numTheories = (correctForContext) ? 4 : 1;
        
        counts = new double[this.numTheories][4];

        sums = new DoubleMatrix1D[this.numTheories][4];
        unscaledCovarianceSums = new DoubleMatrix2D[this.numTheories][4];

        means = new DoubleMatrix1D[this.numTheories][4];
        inverseCovariances = new DoubleMatrix2D[this.numTheories][4];
        norms = new double[this.numTheories][4];

        for (int basePrevIndex = 0; basePrevIndex < this.numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                sums[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                unscaledCovarianceSums[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);

                means[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                inverseCovariances[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);
            }
        }

        alg = new Algebra();
    }


    /**
     * Add a single training point to the model to estimate the means.
     *
     * @param probMatrix     the matrix of probabilities for the base
     * @param fourIntensity  the four raw intensities for the base
     */
    public void addMeanPoint(double[][] probMatrix, double[] fourIntensity) {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double weight = probMatrix[basePrevIndex][baseCurIndex];                

                DoubleMatrix1D weightedChannelIntensities = (DoubleFactory1D.dense).make(fourIntensity);
                weightedChannelIntensities.assign(F.mult(weight));

                sums[basePrevIndex][baseCurIndex].assign(weightedChannelIntensities, F.plus);
                counts[basePrevIndex][baseCurIndex] += weight;
            }
        }

        readyToCall = false;
    }

    /**
     * Add a single training point to the model to estimate the covariances.
     *
     * @param probMatrix     the matrix of probabilities for the base
     * @param fourIntensity  the four raw intensities for the base
     */
    public void addCovariancePoint(double[][] probMatrix, double[] fourIntensity) {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double weight = probMatrix[basePrevIndex][baseCurIndex];

                DoubleMatrix1D mean = sums[basePrevIndex][baseCurIndex].copy();
                mean.assign(F.div(counts[basePrevIndex][baseCurIndex]));

                DoubleMatrix1D sub = (DoubleFactory1D.dense).make(fourIntensity);
                sub.assign(mean, F.minus);

                DoubleMatrix2D cov = (DoubleFactory2D.dense).make(4, 4);
                alg.multOuter(sub, sub, cov);

                cov.assign(F.mult(weight));
                unscaledCovarianceSums[basePrevIndex][baseCurIndex].assign(cov, F.plus);
            }
        }
        
        readyToCall = false;
    }

    /**
     * Precompute all the matrix inversions and determinants we'll need for computing the likelihood distributions.
     */
    public void prepareToCallBases() {
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                means[basePrevIndex][baseCurIndex] = sums[basePrevIndex][baseCurIndex].copy();
                means[basePrevIndex][baseCurIndex].assign(F.div(counts[basePrevIndex][baseCurIndex]));

                inverseCovariances[basePrevIndex][baseCurIndex] = unscaledCovarianceSums[basePrevIndex][baseCurIndex].copy();
                inverseCovariances[basePrevIndex][baseCurIndex].assign(F.div(counts[basePrevIndex][baseCurIndex]));
                DoubleMatrix2D invcov = alg.inverse(inverseCovariances[basePrevIndex][baseCurIndex]);
                inverseCovariances[basePrevIndex][baseCurIndex] = invcov;

                norms[basePrevIndex][baseCurIndex] = Math.pow(alg.det(invcov), 0.5)/Math.pow(2.0*Math.PI, 2.0);
            }
        }

        readyToCall = true;
    }

    /**
     * Compute the likelihood matrix for a base.
     *
     * @param cycle         the cycle we're calling right now
     * @param fourintensity the four intensities of the current cycle's base
     * @return              a 4x4 matrix of likelihoods, where the row is the previous cycle base hypothesis and
     *                      the column is the current cycle base hypothesis
     */
    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        if (!readyToCall) {
            prepareToCallBases();
        }

        double[][] likedist = new double[numTheories][4];
        for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double norm = norms[basePrevIndex][baseCurIndex];

                DoubleMatrix1D sub = (DoubleFactory1D.dense).make(fourintensity);
                sub.assign(means[basePrevIndex][baseCurIndex], F.minus);

                DoubleMatrix1D Ax = alg.mult(inverseCovariances[basePrevIndex][baseCurIndex], sub);
                double exparg = -0.5*alg.mult(sub, Ax);

                likedist[basePrevIndex][baseCurIndex] = norm*Math.exp(exparg);
            }
        }

        return likedist;
    }

    /**
     * Write the model parameters to disk.
     *
     * @param outparam  the file in which the output parameters should be stored
     */
    public void write(File outparam) {
        try {
            PrintWriter writer = new PrintWriter(outparam);

            for (int basePrevIndex = 0; basePrevIndex < numTheories; basePrevIndex++) {
                for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                    writer.print("mean_" + BaseUtils.baseIndexToSimpleBase(baseCurIndex) + " = c(");
                    for (int channel = 0; channel < 4; channel++) {
                        writer.print(sums[basePrevIndex][baseCurIndex].getQuick(channel)/counts[basePrevIndex][baseCurIndex]);

                        if (channel < 3) {
                            writer.print(", ");
                        }
                    }
                    writer.println(");");

                    DoubleMatrix2D cov = unscaledCovarianceSums[basePrevIndex][baseCurIndex].copy();
                    cov.assign(F.div(counts[basePrevIndex][baseCurIndex]));

                    writer.print("cov_" + BaseUtils.baseIndexToSimpleBase(baseCurIndex) + " = matrix(c(");
                    for (int channel1 = 0; channel1 < 4; channel1++) {
                        for (int channel2 = 0; channel2 < 4; channel2++) {
                            writer.print(cov.get(channel2, channel1) + (channel1 == 3 && channel2 == 3 ? "" : ","));
                        }
                    }
                    writer.println("), nr=4, nc=4);\n");
                }
            }

            writer.close();
        } catch (IOException e) {
        }
    }
}
