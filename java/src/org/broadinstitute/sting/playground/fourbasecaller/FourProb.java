package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.QualityUtils;

/**
 * FourProb represents four base hypotheses, their probabilities, and their ranking among one another.
 *
 * @author Kiran Garimella
 */
public class FourProb {
    private double[] baseProbs;
    private int[] baseIndices;

    /**
     * Constructor for FourProb.
     *
     * @param baseIndices the unsorted base indices (A:0, C:1, G:2, T:3).  Now that I think about it, this is stupid.
     * @param baseProbs   the unsorted base hypothesis probabilities.
     */
    public FourProb(int[] baseIndices, double[] baseProbs) {
        Integer[] perm = Utils.SortPermutation(baseProbs);
        double[] ascendingBaseProbs = Utils.PermuteArray(baseProbs, perm);
        int[] ascendingBaseIndices = Utils.PermuteArray(baseIndices, perm);

        this.baseProbs = new double[4];
        this.baseIndices = new int[4];

        for (int i = 0; i < 4; i++) {
            this.baseProbs[i] = ascendingBaseProbs[3 - i];
            this.baseIndices[i] = ascendingBaseIndices[3 - i];
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
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the base label (A, C, G, T).
     */
    public char baseAtRank(int rank) { return baseIndexToBase(indexAtRank(rank)); }

    /**
     * Returns the probability of the base at the specified rank.
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the probability of the base (0.0-1.0)
     */
    public double probAtRank(int rank) { return baseProbs[rank]; }

    /**
     * Returns the quality score of the base at the specified rank.
     * @param rank (0 = best, 3 = worst) the rank of the base whose index should be returned
     * @return the quality score of the base (0-40)
     */
    public byte qualAtRank(int rank) { return QualityUtils.probToQual(probAtRank(rank)); }

    /**
     * A utility method to convert a base index into a base label.
     * @param baseIndex the index of the base (0, 1, 2, 3).
     * @return A, C, G, T, or '.' if the base index can't be understood.
     */
    private char baseIndexToBase(int baseIndex) {
        switch (baseIndex) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '.';
        }
    }

    /**
     * Prettily formats the FourProb info.
     * 
     * @return a prettily formatted Sting containing the base and quality score in rank order.
     */
    public String toString() {
        return (
                "[" + baseAtRank(0) + ":" + qualAtRank(0) + " "
                    + baseAtRank(1) + ":" + qualAtRank(1) + " "
                    + baseAtRank(2) + ":" + qualAtRank(2) + " "
                    + baseAtRank(3) + ":" + qualAtRank(3) + "]"
               );
    }
}
