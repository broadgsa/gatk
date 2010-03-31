package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.Pair;

import java.io.PrintStream;


/**
 * @author aaron
 *
 * Class ValidateDbSNPConversion
 *
 * a quick walker to validate the dbSNP conversion.
 */
public class ValidateDbSNPConversion extends RefWalker<Pair<Matrix.BASE, Matrix.BASE>, Matrix> {
    @Override
    public Pair<Matrix.BASE, Matrix.BASE> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (!tracker.hasROD("dbsnp")) return null;
        rodDbSNP rod =  tracker.lookup("dbSNP",rodDbSNP.class);
        if (rod != null && rod.isSNP() && rod.isBiallelic()) {
            return new Pair<Matrix.BASE, Matrix.BASE>(Matrix.BASE.toBase((byte) ref.getBase()), Matrix.BASE.toBase((byte) rod.getAlternativeBaseForSNP()));
        }
        return null;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public Matrix reduceInit() {
        return new Matrix();
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     *
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public Matrix reduce(Pair<Matrix.BASE, Matrix.BASE> value, Matrix sum) {
        if (value != null) sum.addValue(value.first,value.second);
        return sum;
    }

    public void onTraversalDone(Matrix result) {
        result.printStats(out);
    }

}


class Matrix {
    public enum BASE {
        A((byte) 65), C((byte) 67), G((byte) 71), T((byte) 84), N((byte)78);
        private byte val = 65;

        BASE(byte val) {
            this.val = val;
        }

        byte getVal() {
            return val;
        }

        static BASE toBase(byte b) {
            if (b > 96) b =- 32;
            if (b == A.val) return A;
            if (b == C.val) return C;
            if (b == G.val) return G;
            if (b == T.val) return T;
            // we don't really care, let return N instead of: throw new UnsupportedOperationException("Unknown base " + b);
            return N;
        }
    }

    private int[][] mat = new int[BASE.values().length][BASE.values().length];

    public void addValue(BASE ref, BASE alt) {
        mat[ref.ordinal()][alt.ordinal()] += 1;
    }

    public void printStats(PrintStream stream) {
        for (BASE a : BASE.values()) {
            stream.print("ref " + a.toString() + ":");
            for (BASE b : BASE.values()) {
                stream.print(" " + mat[a.ordinal()][b.ordinal()]);
            }
            stream.println();
        }
    }

    public void addMatrix(Matrix mt) {
        for (BASE a : BASE.values()) {
            for (BASE b : BASE.values()) {
                mat[a.ordinal()][b.ordinal()] += mt.mat[a.ordinal()][b.ordinal()];
            }
        }
    }

    public Matrix() {
        for (int x = 0; x < BASE.values().length; x++) {
            for (int y = 0; y < BASE.values().length; y++) {
                mat[x][y] = 0;
            }
        }
    }
}