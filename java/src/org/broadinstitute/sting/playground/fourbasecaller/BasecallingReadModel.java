package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.Utils;

public class BasecallingReadModel {
    private BasecallingBaseModel[] basemodels = null;

    public BasecallingReadModel(int readLength) {
        basemodels = new BasecallingBaseModel[readLength];

        for (int i = 0; i < readLength; i++) {
            basemodels[i] = new BasecallingBaseModel();
        }
    }

    public void addTrainingPoint(int cycle, char basePrev, char baseCur, byte qualCur, double[] fourintensity) {
        basemodels[cycle].addTrainingPoint(basePrev, baseCur, qualCur, fourintensity);
    }

    public double[][] computeLikelihoods(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, basePrev, qualPrev, fourintensity);
    }

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
