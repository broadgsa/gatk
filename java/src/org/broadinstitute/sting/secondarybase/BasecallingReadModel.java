package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;
import java.util.ArrayList;

/**
 * BasecallingReadModel represents the statistical models for
 * all bases in all cycles.  It allows for easy training via
 * the addTrainingPoint() method, and for the computation of
 * the 4x4 likelihood matrix or the 1x4 probability vector
 * (with contextual components marginalized out of the
 * likelihood matrix).
 *
 * @author Kiran Garimella
 */
public class BasecallingReadModel {
    private BasecallingBaseModel[] basemodels = null;
    private boolean correctForContext = true;

    public BasecallingReadModel(int readLength) {
        initialize(readLength);
    }

    public BasecallingReadModel(ArrayList<RawRead> trainingData) {
        initialize(trainingData.get(0).getReadLength());

        train(trainingData);
    }

    public void initialize(int readLength) {
        basemodels = new BasecallingBaseModel[readLength];

        for (int cycle = 0; cycle < readLength; cycle++) {
            basemodels[cycle] = new BasecallingBaseModel(cycle != 0 && correctForContext);
        }
    }

    public void train(ArrayList<RawRead> trainingData) {
        //for (int readIndex = 0; readIndex < trainingData.size(); readIndex++) {
        for ( RawRead read : trainingData ) {
            addMeanPoints(read);
        }

        for ( RawRead read : trainingData ) {
            addCovariancePoints(read);
        }
    }

    public void addMeanPoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addMeanPoint(probMatrix, fourintensity);
    }

    public void addMeanPoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) ints[cycle][channel];
            }

            basemodels[cycle].addMeanPoint(probMatrix, fourIntensity);
        }
    }

    public void addCovariancePoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addCovariancePoint(probMatrix, fourintensity);
    }

    public void addCovariancePoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) ints[cycle][channel];
            }

            basemodels[cycle].addCovariancePoint(probMatrix, fourIntensity);
        }
    }

    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, fourintensity);
    }

    public FourProb computeProbabilities(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        double[][] likes = computeLikelihoods(cycle, fourintensity);

        double total = 0;

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double prior = 1.0;
                if (correctForContext) {
                    double prob = QualityUtils.qualToProb(qualPrev);
                    if (basePrevIndex == BaseUtils.simpleBaseToBaseIndex(basePrev)) {
                        prior = prob;
                    } else {
                        prior = (1.0 - prob)/((double) (4*likes.length - 1));
                    }
                }
                likes[basePrevIndex][baseCurIndex] = prior*likes[basePrevIndex][baseCurIndex];
                total += likes[basePrevIndex][baseCurIndex];
            }
        }

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                likes[basePrevIndex][baseCurIndex] /= total;
            }
        }

        return new FourProb(likes);
    }

    public FourProbRead call(RawRead read) {
        FourProbRead fpr = new FourProbRead(read.getReadLength());
        
        for (int cycle = 0; cycle < read.getReadLength(); cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : read.getSequence()[cycle - 1]);
            byte qualPrev = ((cycle == 0) ? 0 : read.getQuals()[cycle - 1]);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) read.getIntensities()[cycle][channel];
            }

            fpr.add(cycle, computeProbabilities(cycle, basePrev, qualPrev, fourIntensity));

        }

        return fpr;
    }

    public double[][] getBaseProbabilityMatrix(int cycle, char basePrev, char baseCur, double probCur) {
        double[][] dist = new double[(correctForContext && cycle > 0) ? 4 : 1][4];

        int actualBasePrevIndex = (correctForContext && cycle > 0) ? BaseUtils.simpleBaseToBaseIndex(basePrev) : 0;
        int actualBaseCurIndex = BaseUtils.simpleBaseToBaseIndex(baseCur);

        double residualTheories = (double) (dist.length*dist[0].length - 1);

        for (int basePrevIndex = 0; basePrevIndex < dist.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < dist[basePrevIndex].length; baseCurIndex++) {
                dist[basePrevIndex][baseCurIndex] = (basePrevIndex == actualBasePrevIndex && baseCurIndex == actualBaseCurIndex) ? probCur : ((1.0 - probCur)/residualTheories);
            }
        }

        return dist;
    }

    public void write(File dir) {
        for (int cycle = 0; cycle < basemodels.length; cycle++) {
            File outparam = new File(dir.getPath() + "/param." + cycle + ".r");
            basemodels[cycle].write(outparam);
        }
    }
}
