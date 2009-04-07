package org.broadinstitute.sting.playground.fourbasecaller;

import java.io.File;

/**
 * BasecallingReadModel represents the statistical models for
 * all bases in all cycles.  It allows for easy, one-pass
 * training via the addTrainingPoint() method, and for the
 * computation of the 4x4 likelihood matrix or the 1x4
 * probability vector (with contextual components marginalized
 * out of the likelihood matrix).
 *
 * @author Kiran Garimella
 */
public class BasecallingReadModel {
    private BasecallingBaseModel[] basemodels = null;

    /**
     * Constructor for BasecallingReadModel.
     *
     * @param readLength    the length of the read that this model will support
     */
    public BasecallingReadModel(int readLength) {
        basemodels = new BasecallingBaseModel[readLength];

        for (int i = 0; i < readLength; i++) {
            basemodels[i] = new BasecallingBaseModel();
        }
    }

    /**
     * Add a single training point to the model.
     *
     * @param cycle         the cycle for which this point should be added
     * @param baseCur       the current base
     * @param qualCur       the current base's quality
     * @param fourintensity the four intensities of the current base
     */
    public void addMeanPoint(int cycle, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addMeanPoint(baseCur, qualCur, fourintensity);
    }

    public void addCovariancePoint(int cycle, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addCovariancePoint(baseCur, qualCur, fourintensity);
    }

    /**
     * Compute the likelihood matrix for a given cycle.
     *
     * @param cycle         the cycle number for the current base
     * @param fourintensity the four intensities for the current cycle's base
     * @return 4x4 matrix of likelihoods
     */
    public double[] computeLikelihoods(int cycle, double[] fourintensity) {
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
    public FourProb computeProbabilities(int cycle, double[] fourintensity) {
        double[] likes = computeLikelihoods(cycle, fourintensity);

        double total = 0;

        for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) { total += likes[baseCurIndex]; }
        for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) { likes[baseCurIndex] /= total; }

        return new FourProb(likes);
    }

    /**
     * Writes model parameters to a file per cycle.
     *
     * @param dir  the directory where the parameters should be written
     */
    public void write(File dir) {
        for (int cycle = 0; cycle < basemodels.length; cycle++) {
            File outparam = new File(dir.getPath() + "/param." + cycle);
            basemodels[cycle].write(outparam);
        }
    }
}
