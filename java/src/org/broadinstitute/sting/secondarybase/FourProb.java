package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;

/**
 * FourProb represents four base hypotheses, their probabilities, and the ranking among one another.
 *
 * @author Kiran Garimella
 */
public class FourProb {
    private int[] baseIndices;
    private double[] baseProbs;

    /**
     * Constructor for FourProb.
     *
     * @param baseLikes  the unsorted base hypothesis probabilities (in ACGT order).
     */
    public FourProb(double[][] baseLikes) {
        double[] baseProbs = new double[4];
        for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
            for (int basePrevIndex = 0; basePrevIndex < baseLikes.length; basePrevIndex++) {
                baseProbs[baseCurIndex] += baseLikes[basePrevIndex][baseCurIndex];
            }
        }

        this.baseProbs = new double[4];
        this.baseIndices = new int[4];

        // store this information in sorted form
        for (int i = 0; i < 4; i++) {
            if (baseProbs[i] > this.baseProbs[0]) {
                this.baseProbs[1] = this.baseProbs[0];
                this.baseIndices[1] = this.baseIndices[0];

                this.baseProbs[0] = baseProbs[i];
                this.baseIndices[0] = i;
            } else if (baseProbs[i] > this.baseProbs[1]) {
                this.baseProbs[2] = this.baseProbs[1];
                this.baseIndices[2] = this.baseIndices[1];

                this.baseProbs[1] = baseProbs[i];
                this.baseIndices[1] = i;
            } else if (baseProbs[i] > this.baseProbs[2]) {
                this.baseProbs[3] = this.baseProbs[2];
                this.baseIndices[3] = this.baseIndices[2];

                this.baseProbs[2] = baseProbs[i];
                this.baseIndices[2] = i;
            } else {
                this.baseProbs[3] = baseProbs[i];
                this.baseIndices[3] = i;
            }
        }
    }

    /**
     * Returns the index of the base at the specified rank.
     *
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the index (0, 1, 2, 3).
     */
    public int indexAtRank(int rank) { return baseIndices[rank]; }

    /**
     * Returns the base label of the base at the specified rank.
     *
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the base label (A, C, G, T).
     */
    //public char baseAtRank(int rank) { return baseIndexToBase(indexAtRank(rank)); }
    public char baseAtRank(int rank) { return BaseUtils.baseIndexToSimpleBase(indexAtRank(rank)); }

    /**
     * Returns the probability of the base at the specified rank.
     *
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the probability of the base (0.0-1.0)
     */
    public double probAtRank(int rank) { return baseProbs[rank]; }

    /**
     * Returns the quality score of the base at the specified rank.
     *
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the quality score of the base (0-40)
     */
    public byte qualAtRank(int rank) { return QualityUtils.probToQual(probAtRank(rank)); }

    /**
     * Prettily formats the FourProb info.
     * 
     * @return a prettily formatted Sting containing the base and quality score in rank order.
     */
    public String toString() {
        return (
                "[" + baseAtRank(0) + ":" + probAtRank(0) + " "
                    + baseAtRank(1) + ":" + probAtRank(1) + " "
                    + baseAtRank(2) + ":" + probAtRank(2) + " "
                    + baseAtRank(3) + ":" + probAtRank(3) + "]"
               );
    }
}
