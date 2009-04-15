package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

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

    /**
     * Add a single training point to the model means.
     *
     * @param cycle         the cycle for which this point should be added
     * @param basePrev      the previous base
     * @param baseCur       the current base
     * @param qualCur       the current base's quality
     * @param fourintensity the four intensities of the current base
     */
    public void addMeanPoint(int cycle, char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addMeanPoint(basePrev, baseCur, qualCur, fourintensity);
    }

    /**
     * Add a single training point to the model covariances.
     *
     * @param cycle         the cycle for which this point should be added
     * @param basePrev      the previous base
     * @param baseCur       the current base
     * @param qualCur       the current base's quality
     * @param fourintensity the four intensities of the current base
     */
    public void addCovariancePoint(int cycle, char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addCovariancePoint(basePrev, baseCur, qualCur, fourintensity);
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

        /*
        System.out.println(fp.baseAtRank(0) + ":");
        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                System.out.print(likes[basePrevIndex][baseCurIndex] + " ");
            }
            System.out.println();
        }
        System.out.println();
        */

        //return new FourProb(likes);
        return fp;
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
