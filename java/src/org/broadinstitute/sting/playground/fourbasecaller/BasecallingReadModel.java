package org.broadinstitute.sting.playground.fourbasecaller;

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
     * @param cycle         the cycle for which this point should be added
     * @param basePrev      the previous base
     * @param baseCur       the current base
     * @param qualCur       the current base's quality
     * @param fourintensity the four intensities of the current base
     */
    public void addTrainingPoint(int cycle, char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addTrainingPoint(basePrev, baseCur, qualCur, fourintensity);
    }

    /**
     * Compute the likelihood matrix for a given cycle.
     *
     * @param cycle         the cycle number for the current base
     * @param basePrev      the previous cycle's base
     * @param qualPrev      the quality score for the previous cycle's base
     * @param fourintensity the four intensities for the current cycle's base
     * @return 4x4 matrix of likelihoods
     */
    public double[][] computeLikelihoods(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, basePrev, qualPrev, fourintensity);
    }

    /**
     * Compute the probability distribution for the base at a given cycle.
     * Contextual components of the likelihood matrix are marginalized out.
     *
     * @param cycle         the cycle number for the current base
     * @param basePrev      the previous cycle's base
     * @param qualPrev      the quality score for the previous cycle's base
     * @param fourintensity the four intensities for the current cycle's base
     * @return an instance of FourProb, which encodes a base hypothesis, its probability,
     *         and the ranking among the other hypotheses
     */
    public FourProb computeProbabilities(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        double[][] likes = computeLikelihoods(cycle, basePrev, qualPrev, fourintensity);

        double[] probs = new double[4];
        int[] baseindices = { 0, 1, 2, 3 };
        double total = 0;

        for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
            for (int basePrevIndex = 0; basePrevIndex < 4; basePrevIndex++) {
                probs[baseCurIndex] += likes[basePrevIndex][baseCurIndex];
            }
            total += probs[baseCurIndex];
        }

        for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
            probs[baseCurIndex] /= total;
        }

        return new FourProb(baseindices, probs);
    }
}
