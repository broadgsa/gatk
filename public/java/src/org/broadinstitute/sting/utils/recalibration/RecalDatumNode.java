package org.broadinstitute.sting.utils.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.HashSet;
import java.util.Set;

/**
 * A tree of recal datum, where each contains a set of sub datum representing sub-states of the higher level one
 *
 * @author Mark DePristo
 * @since 07/27/12
 */
public class RecalDatumNode<T extends RecalDatum> {
    protected static Logger logger = Logger.getLogger(RecalDatumNode.class);
    private final static double UNINITIALIZED = -1.0;
    private final T recalDatum;
    private double fixedPenalty = UNINITIALIZED;
    private final Set<RecalDatumNode<T>> subnodes;

    public RecalDatumNode(final T recalDatum) {
        this(recalDatum, new HashSet<RecalDatumNode<T>>());
    }

    @Override
    public String toString() {
        return recalDatum.toString();
    }

    public RecalDatumNode(final T recalDatum, final Set<RecalDatumNode<T>> subnodes) {
        this(recalDatum, UNINITIALIZED, subnodes);
    }

    protected RecalDatumNode(final T recalDatum, final double fixedPenalty) {
        this(recalDatum, fixedPenalty, new HashSet<RecalDatumNode<T>>());
    }

    protected RecalDatumNode(final T recalDatum, final double fixedPenalty, final Set<RecalDatumNode<T>> subnodes) {
        this.recalDatum = recalDatum;
        this.fixedPenalty = fixedPenalty;
        this.subnodes = new HashSet<RecalDatumNode<T>>(subnodes);
    }

    public T getRecalDatum() {
        return recalDatum;
    }

    public Set<RecalDatumNode<T>> getSubnodes() {
        return subnodes;
    }

    public double getPenalty() {
        if ( fixedPenalty != UNINITIALIZED )
            return fixedPenalty;
        else
            return calcPenalty(recalDatum.getEmpiricalErrorRate());
    }

    public double calcAndSetFixedPenalty(final boolean doEntireTree) {
        fixedPenalty = calcPenalty(recalDatum.getEmpiricalErrorRate());
        if ( doEntireTree )
            for ( final RecalDatumNode<T> sub : subnodes )
                sub.calcAndSetFixedPenalty(doEntireTree);
        return fixedPenalty;
    }

    public void addSubnode(final RecalDatumNode<T> sub) {
        subnodes.add(sub);
    }

    public boolean isLeaf() {
        return subnodes.isEmpty();
    }

    public int getNumBranches() {
        return subnodes.size();
    }

    public double getMinNodePenalty() {
        if ( isLeaf() )
            return Double.MAX_VALUE;
        else {
            double minPenalty = getPenalty();
            for ( final RecalDatumNode<T> sub : subnodes )
                minPenalty = Math.min(minPenalty, sub.getMinNodePenalty());
            return minPenalty;
        }
    }

    public int maxDepth() {
        int subMax = 0;
        for ( final RecalDatumNode<T> sub : subnodes )
            subMax = Math.max(subMax, sub.maxDepth());
        return subMax + 1;
    }

    public int size() {
        int size = 1;
        for ( final RecalDatumNode<T> sub : subnodes )
            size += sub.size();
        return size;
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
            return (Math.abs(Math.log10(recalDatum.getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * recalDatum.getNumObservations();
            // TODO -- how we can generalize this calculation?
//            if ( this.qEnd <= minInterestingQual )
//                // It's free to merge up quality scores below the smallest interesting one
//                return 0;
//            else {
//                return (Math.abs(Math.log10(getEmpiricalErrorRate()) - Math.log10(globalErrorRate))) * getNumObservations();
//            }
        } else {
            double sum = 0;
            for ( final RecalDatumNode<T> hrd : subnodes)
                sum += hrd.calcPenalty(globalErrorRate);
            return sum;
        }
    }

    public RecalDatumNode<T> pruneToDepth(final int maxDepth) {
        if ( maxDepth < 1 )
            throw new IllegalArgumentException("maxDepth < 1");
        else {
            final Set<RecalDatumNode<T>> subPruned = new HashSet<RecalDatumNode<T>>(getNumBranches());
            if ( maxDepth > 1 )
                for ( final RecalDatumNode<T> sub : subnodes )
                    subPruned.add(sub.pruneToDepth(maxDepth - 1));
            return new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, subPruned);
        }
    }

    public RecalDatumNode<T> pruneByPenalty(final int maxElements) {
        RecalDatumNode<T> root = this;

        while ( root.size() > maxElements ) {
            // remove the lowest penalty element, and continue
            root = root.removeLowestPenaltyNode();
        }

        // our size is below the target, so we are good, return
        return root;
    }

    /**
     * Find the lowest penalty node in the tree, and return a tree without it
     *
     * Note this excludes the current (root) node
     *
     * @return
     */
    private RecalDatumNode<T> removeLowestPenaltyNode() {
        final RecalDatumNode<T> oneRemoved = removeFirstNodeWithPenalty(getMinNodePenalty()).getFirst();
        if ( oneRemoved == null )
            throw new IllegalStateException("Removed our root node, wow, didn't expect that");
        return oneRemoved;
    }

    private Pair<RecalDatumNode<T>, Boolean> removeFirstNodeWithPenalty(final double penaltyToRemove) {
        if ( getPenalty() == penaltyToRemove ) {
            logger.info("Removing " + this + " with penalty " + penaltyToRemove);
            if ( isLeaf() )
                throw new IllegalStateException("Trying to remove a leaf node from the tree! " + this + " " + penaltyToRemove);
            // node is the thing we are going to remove, but without any subnodes
            final RecalDatumNode<T> node = new RecalDatumNode<T>(getRecalDatum(), fixedPenalty);
            return new Pair<RecalDatumNode<T>, Boolean>(node, true);
        } else {
            // did we remove something in a sub branch?
            boolean removedSomething = false;

            // our sub nodes with the penalty node removed
            final Set<RecalDatumNode<T>> sub = new HashSet<RecalDatumNode<T>>(getNumBranches());

            for ( final RecalDatumNode<T> sub1 : subnodes ) {
                if ( removedSomething ) {
                    // already removed something, just add sub1 back to sub
                    sub.add(sub1);
                } else {
                    // haven't removed anything yet, so try
                    final Pair<RecalDatumNode<T>, Boolean> maybeRemoved = sub1.removeFirstNodeWithPenalty(penaltyToRemove);
                    removedSomething = maybeRemoved.getSecond();
                    sub.add(maybeRemoved.getFirst());
                }
            }

            final RecalDatumNode<T> node = new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, sub);
            return new Pair<RecalDatumNode<T>, Boolean>(node, removedSomething);
        }
    }
}
