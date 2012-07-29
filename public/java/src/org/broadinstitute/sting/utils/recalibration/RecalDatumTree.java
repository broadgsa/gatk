package org.broadinstitute.sting.utils.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * A tree of recal datum, where each contains a set of sub datum representing sub-states of the higher level one
 *
 * @author Mark DePristo
 * @since 07/27/12
 */
public class RecalDatumTree extends RecalDatum {
    final Set<RecalDatumTree> subnodes;

    protected RecalDatumTree(final long nObservations, final long nErrors, final byte reportedQual) {
        this(nObservations, nErrors, reportedQual, new HashSet<RecalDatumTree>());
    }

    public RecalDatumTree(final long nObservations, final long nErrors, final byte reportedQual, final Set<RecalDatumTree> subnodes) {
        super(nObservations, nErrors, reportedQual);
        this.subnodes = new HashSet<RecalDatumTree>(subnodes);
    }

    public double getPenalty() {
        return calcPenalty(getEmpiricalErrorRate());
    }

    public void addSubnode(final RecalDatumTree sub) {
        subnodes.add(sub);
    }

    public boolean isLeaf() {
        return subnodes.isEmpty();
    }

    /**
     * Calculate the penalty of this interval, given the overall error rate for the interval
     *
     * If the globalErrorRate is e, this value is:
     *
     * sum_i |log10(e_i) - log10(e)| * nObservations_i
     *
     * each the index i applies to all leaves of the tree accessible from this interval
     * (found recursively from subnodes as necessary)
     *
     * @param globalErrorRate overall error rate in real space against which we calculate the penalty
     * @return the cost of approximating the bins in this interval with the globalErrorRate
     */
    @Requires("globalErrorRate >= 0.0")
    @Ensures("result >= 0.0")
    private double calcPenalty(final double globalErrorRate) {
        if ( globalErrorRate == 0.0 ) // there were no observations, so there's no penalty
            return 0.0;

        if ( isLeaf() ) {
            // this is leave node
            return (Math.abs(Math.log10(getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * getNumObservations();
            // TODO -- how we can generalize this calculation?
//            if ( this.qEnd <= minInterestingQual )
//                // It's free to merge up quality scores below the smallest interesting one
//                return 0;
//            else {
//                return (Math.abs(Math.log10(getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * getNumObservations();
//            }
        } else {
            double sum = 0;
            for ( final RecalDatumTree hrd : subnodes)
                sum += hrd.calcPenalty(globalErrorRate);
            return sum;
        }
    }
}
