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
    private boolean correctForContext = false;

    /**
     * Constructor for BasecallingReadModel.
     *
     * @param readLength    the length of the read that this model will support
     */
    public BasecallingReadModel(int readLength, boolean correctForContext) {
        this.correctForContext = correctForContext;
        
        basemodels = new BasecallingBaseModel[readLength];

        for (int cycle = 0; cycle < readLength; cycle++) {
            basemodels[cycle] = new BasecallingBaseModel(cycle == 0 ? false : correctForContext);
        }
    }

    public void train(BasecallingTrainingSet trainingSet) {
        ArrayList<RawRead> trainingData = trainingSet.getTrainingData();

        for (int readIndex = 0; readIndex < trainingData.size(); readIndex++) {
            RawRead read = trainingData.get(readIndex);
            addMeanPoints(read);
        }

        for (int readIndex = 0; readIndex < trainingData.size(); readIndex++) {
            RawRead read = trainingData.get(readIndex);
            addCovariancePoints(read);
        }
    }

    /**
     * Add a single training point to the model means.
     *
     * @param cycle         the cycle for which this point should be added
     * @param fourintensity the four intensities of the current base
     */
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

    /**
     * Add a single training point to the model covariances.
     *
     * @param cycle         the cycle for which this point should be added
     * @param fourintensity the four intensities of the current base
     */
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

    /**
     * Compute the likelihood matrix for a given cycle.
     *
     * @param cycle         the cycle number for the current base
     * @param fourintensity the four intensities for the current cycle's base
     * @return 4x4 matrix of likelihoods
     */
    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, fourintensity);
    }

    /**
     * Compute the probability distribution for the base at a given cycle.
     *
     * @param cycle         the cycle number for the current base
     * @param fourintensity the four intensities for the current cycle's base
     * @return an instance of FourProb, which encodes a base hypothesis, its probability,
     *         and the ranking among the other hypotheses
     */
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

        FourProb fp = new FourProb(likes);

        return fp;
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

            //fps[cycle] = computeProbabilities(cycle, basePrev, qualPrev, fourIntensity);
            fpr.add(cycle, computeProbabilities(cycle, basePrev, qualPrev, fourIntensity));

        }

        return fpr;
    }

    /**
     * Returns the base probability matrix
     *
     * @param cycle        the cycle of the base
     * @param basePrev     the previous base
     * @param baseCur      the current base
     * @param probCur      the probability of the current base
     * @return the base probability matrix (1x4 if no contextual correction, 4x4 otherwise)
     */
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

    /**
     * Writes model parameters to a file per cycle.
     *
     * @param dir  the directory where the parameters should be written
     */
    public void write(File dir) {
        for (int cycle = 0; cycle < basemodels.length; cycle++) {
            File outparam = new File(dir.getPath() + "/param." + cycle + ".r");
            basemodels[cycle].write(outparam);
        }
    }
}
