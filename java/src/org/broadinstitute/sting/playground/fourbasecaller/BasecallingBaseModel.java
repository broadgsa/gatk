package org.broadinstitute.sting.playground.fourbasecaller;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.linalg.Algebra;

import org.broadinstitute.sting.utils.QualityUtils;

import java.io.*;

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

    private DoubleMatrix1D[][] runningChannelSums;
    private DoubleMatrix2D[][] runningChannelProductSums;

    private boolean readyToCall = false;
    private DoubleMatrix1D[][] means;
    private DoubleMatrix2D[][] inverseCovariances;
    private double[][] norms;

    private Algebra alg;

    /**
     * Constructor for BasecallingBaseModel
     */
    public BasecallingBaseModel() {
        counts = new double[4][4];

        runningChannelSums = new DoubleMatrix1D[4][4];
        runningChannelProductSums = new DoubleMatrix2D[4][4];

        means = new DoubleMatrix1D[4][4];
        inverseCovariances = new DoubleMatrix2D[4][4];
        norms = new double[4][4];

        for (int basePrevIndex = 0; basePrevIndex < 4; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                runningChannelSums[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                runningChannelProductSums[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);

                means[basePrevIndex][baseCurIndex] = (DoubleFactory1D.dense).make(4);
                inverseCovariances[basePrevIndex][baseCurIndex] = (DoubleFactory2D.dense).make(4, 4);
            }
        }

        alg = new Algebra();
    }

    /**
     * Add a single training point to the model.
     *
     * @param basePrev       the previous cycle's base call (A, C, G, T, or * for the first cycle)
     * @param baseCur        the current cycle's base call (A, C, G, T)
     * @param qualCur        the quality score for the current cycle's base call
     * @param fourintensity  the four intensities for the current cycle's base call
     */
    public void addTrainingPoint(char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        int actualBasePrevIndex = baseToBaseIndex(basePrev);
        int actualBaseCurIndex = baseToBaseIndex(baseCur);
        double actualWeight = QualityUtils.qualToProb(qualCur);
        double otherTheories = (basePrev == '*') ? 3.0 : 15.0;

        cern.jet.math.Functions F = cern.jet.math.Functions.functions;

        for (int basePrevIndex = 0; basePrevIndex < 4; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                // We want to upweight the correct theory as much as we can and spread the remainder out evenly between all other hypotheses.
                double weight = (basePrevIndex == actualBasePrevIndex && baseCurIndex == actualBaseCurIndex) ? actualWeight : ((1.0 - actualWeight)/otherTheories);

                DoubleMatrix1D weightedChannelIntensities = (DoubleFactory1D.dense).make(fourintensity);
                weightedChannelIntensities.assign(F.mult(weight));

                runningChannelSums[basePrevIndex][baseCurIndex].assign(weightedChannelIntensities, F.plus);
                counts[basePrevIndex][baseCurIndex] += weight;
            }
        }

        /*
        if (basePrevIndex >= 0 && baseCurIndex >= 0) {
            for (int channel = 0; channel < 4; channel++) {
                double weight = QualityUtils.qualToProb(qualCur);
                double channelIntensity = fourintensity[channel];

                runningChannelSums[basePrevIndex][baseCurIndex].setQuick(channel, runningChannelSums[basePrevIndex][baseCurIndex].getQuick(channel) + weight*channelIntensity);

                for (int cochannel = 0; cochannel < 4; cochannel++) {
                    double cochannelIntensity = fourintensity[cochannel];
                    runningChannelProductSums[basePrevIndex][baseCurIndex].setQuick(channel, cochannel, runningChannelProductSums[basePrevIndex][baseCurIndex].getQuick(channel, cochannel) + weight*channelIntensity*cochannelIntensity);
                }
            }

            counts[basePrevIndex][baseCurIndex]++;
        }
        */

        readyToCall = false;
    }

    /**
     * Precompute all the matrix inversions and determinants we'll need for computing the likelihood distributions.
     */
    public void prepareToCallBases() {
        for (int basePrevIndex = 0; basePrevIndex < 4; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                for (int channel = 0; channel < 4; channel++) {
                    means[basePrevIndex][baseCurIndex].setQuick(channel, runningChannelSums[basePrevIndex][baseCurIndex].getQuick(channel)/counts[basePrevIndex][baseCurIndex]);
                    
                    for (int cochannel = 0; cochannel < 4; cochannel++) {
                        // Cov(Xi, Xj) = E(XiXj) - E(Xi)E(Xj)
                        inverseCovariances[basePrevIndex][baseCurIndex].setQuick(channel, cochannel, (runningChannelProductSums[basePrevIndex][baseCurIndex].getQuick(channel, cochannel)/counts[basePrevIndex][baseCurIndex]) - (runningChannelSums[basePrevIndex][baseCurIndex].getQuick(channel)/counts[basePrevIndex][baseCurIndex])*(runningChannelSums[basePrevIndex][baseCurIndex].getQuick(cochannel)/counts[basePrevIndex][baseCurIndex]));
                    }
                }

                DoubleMatrix2D invcov = alg.inverse(inverseCovariances[basePrevIndex][baseCurIndex]);
                inverseCovariances[basePrevIndex][baseCurIndex] = invcov;
                norms[basePrevIndex][baseCurIndex] = Math.pow(alg.det(invcov), 0.5)/Math.pow(2.0*Math.PI, 2.0);
            }
        }

        readyToCall = true;
    }

    /**
     * Compute the likelihood matrix for a base (contextual priors included).
     *
     * @param cycle         the cycle we're calling right now
     * @param basePrev      the previous cycle's base
     * @param qualPrev      the previous cycle's quality score
     * @param fourintensity the four intensities of the current cycle's base
     * @return              a 4x4 matrix of likelihoods, where the row is the previous cycle base hypothesis and
     *                      the column is the current cycle base hypothesis
     */
    public double[][] computeLikelihoods(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        if (!readyToCall) {
            prepareToCallBases();
        }

        double[][] probdist = new double[4][4];
        double probPrev = (cycle == 0) ? 1.0 : QualityUtils.qualToProb(qualPrev);
        int baseIndex = (cycle == 0) ? 0 : baseToBaseIndex(basePrev);

        for (int basePrevIndex = 0; basePrevIndex < ((cycle == 0) ? 1 : 4); basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double[] diff = new double[4];
                for (int channel = 0; channel < 4; channel++) {
                    diff[channel] = fourintensity[channel] - means[basePrevIndex][baseCurIndex].getQuick(channel);
                }
                
                DoubleMatrix1D sub = (DoubleFactory1D.dense).make(diff);
                DoubleMatrix1D Ax = alg.mult(inverseCovariances[basePrevIndex][baseCurIndex], sub);

                double exparg = -0.5*alg.mult(sub, Ax);
                probdist[basePrevIndex][baseCurIndex] = (baseIndex == basePrevIndex ? probPrev : 1.0 - probPrev)*norms[basePrevIndex][baseCurIndex]*Math.exp(exparg);
            }
        }

        return probdist;
    }

    public void write(File outparam) {
        try {
            PrintWriter writer = new PrintWriter(outparam);

            for (int basePrevIndex = 0; basePrevIndex < 4; basePrevIndex++) {
                for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                    writer.print("mean_" + baseIndexToBase(basePrevIndex) + "" + baseIndexToBase(baseCurIndex) + " : [ ");
                    for (int channel = 0; channel < 4; channel++) {
                        writer.print(runningChannelSums[basePrevIndex][baseCurIndex].getQuick(channel)/counts[basePrevIndex][baseCurIndex]);
                        writer.print(" ");
                    }
                    writer.print("]\n");
                }
            }

            writer.close();
        } catch (IOException e) {

        }
    }

    /**
     * Utility method for converting a base ([Aa*], [Cc], [Gg], [Tt]) to an index (0, 1, 2, 3);
     * @param base
     * @return 0, 1, 2, 3, or -1 if the base can't be understood.
     */
    private int baseToBaseIndex(char base) {
        switch (base) {
            case 'A':
            case 'a':
            case '*': return 0;

            case 'C':
            case 'c': return 1;

            case 'G':
            case 'g': return 2;

            case 'T':
            case 't': return 3;
        }

        return -1;
    }

    private char baseIndexToBase(int baseIndex) {
        switch (baseIndex) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '.';
        }
    }
}
