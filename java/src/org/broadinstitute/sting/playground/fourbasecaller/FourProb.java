package org.broadinstitute.sting.playground.fourbasecaller;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.QualityUtils;

public class FourProb {
    private double[] baseProbs;
    private int[] baseIndices;

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

    public int indexAtRank(int rank) { return baseIndices[rank]; }
    public char baseAtRank(int rank) { return baseIndexToBase(indexAtRank(rank)); }
    public double probAtRank(int rank) { return baseProbs[rank]; }
    public byte qualAtRank(int rank) { return QualityUtils.probToQual(probAtRank(rank)); }

    private char baseIndexToBase(int baseIndex) {
        switch (baseIndex) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '.';
        }
    }

    public String toString() {
        return (
                "[" + baseAtRank(0) + ":" + qualAtRank(0) + " "
                    + baseAtRank(1) + ":" + qualAtRank(1) + " "
                    + baseAtRank(2) + ":" + qualAtRank(2) + " "
                    + baseAtRank(3) + ":" + qualAtRank(3) + "]"
               );
    }
}
