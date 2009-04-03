package org.broadinstitute.sting.playground.fourbasecaller;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.linalg.Algebra;

import org.broadinstitute.sting.utils.QualityUtils;

public class BasecallingBaseModel {
    private double[][] counts;

    private DoubleMatrix1D[][] runningChannelSums;
    private DoubleMatrix2D[][] runningChannelProductSums;

    private boolean readyToCall = false;
    private DoubleMatrix1D[][] means;
    private DoubleMatrix2D[][] inverseCovariances;
    private double[][] norms;

    Algebra alg;

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

    public void addTrainingPoint(char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        int basePrevIndex = baseToBaseIndex(basePrev);
        int baseCurIndex = baseToBaseIndex(baseCur);

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

        readyToCall = false;
    }

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
}
