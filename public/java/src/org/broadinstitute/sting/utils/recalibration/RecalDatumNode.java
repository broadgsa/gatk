package org.broadinstitute.sting.utils.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.MathUtils;
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
    private final static boolean USE_CHI2 = true;
    protected static Logger logger = Logger.getLogger(RecalDatumNode.class);
    private final static double UNINITIALIZED = Double.NEGATIVE_INFINITY;
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
            return calcPenalty();
    }

    public double calcAndSetFixedPenalty(final boolean doEntireTree) {
        fixedPenalty = calcPenalty();
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

    /**
     * Total penalty is the sum of leaf node penalties
     *
     * This algorithm assumes that penalties have been fixed before pruning, as leaf nodes by
     * definition have 0 penalty unless they represent a pruned tree with underlying -- but now
     * pruned -- subtrees
     *
     * @return
     */
    public double totalPenalty() {
        if ( isLeaf() )
            return getPenalty();
        else {
            double sum = 0.0;
            for ( final RecalDatumNode<T> sub : subnodes )
                sum += sub.totalPenalty();
            return sum;
        }
    }

    public int maxDepth() {
        int subMax = 0;
        for ( final RecalDatumNode<T> sub : subnodes )
            subMax = Math.max(subMax, sub.maxDepth());
        return subMax + 1;
    }

    public int minDepth() {
        if ( isLeaf() )
            return 1;
        else {
            int subMin = Integer.MAX_VALUE;
            for ( final RecalDatumNode<T> sub : subnodes )
                subMin = Math.min(subMin, sub.minDepth());
            return subMin + 1;
        }
    }

    public int size() {
        int size = 1;
        for ( final RecalDatumNode<T> sub : subnodes )
            size += sub.size();
        return size;
    }

    public int numLeaves() {
        if ( isLeaf() )
            return 1;
        else {
            int size = 0;
            for ( final RecalDatumNode<T> sub : subnodes )
                size += sub.numLeaves();
            return size;
        }
    }

    private double calcPenalty() {
        if ( USE_CHI2 )
            return calcPenaltyChi2();
        else
            return calcPenaltyLog10(getRecalDatum().getEmpiricalErrorRate());
    }

    private double calcPenaltyChi2() {
        if ( isLeaf() )
            return 0.0;
        else {
            final long[][] counts = new long[subnodes.size()][2];

            int i = 0;
            for ( RecalDatumNode<T> subnode : subnodes ) {
                counts[i][0] = subnode.getRecalDatum().getNumMismatches();
                counts[i][1] = subnode.getRecalDatum().getNumObservations();
                i++;
            }

            final double chi2 = new ChiSquareTestImpl().chiSquare(counts);

//            StringBuilder x = new StringBuilder();
//            StringBuilder y = new StringBuilder();
//            for ( int k = 0; k < counts.length; k++) {
//                if ( k != 0 ) {
//                    x.append(", ");
//                    y.append(", ");
//                }
//                x.append(counts[k][0]);
//                y.append(counts[k][1]);
//            }
//            logger.info("x = c(" + x.toString() + ")");
//            logger.info("y = c(" + y.toString() + ")");
//            logger.info("chi2 = " + chi2);

            return chi2;
            //return Math.log10(chi2);
        }
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
    private double calcPenaltyLog10(final double globalErrorRate) {
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
                sum += hrd.calcPenaltyLog10(globalErrorRate);
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
        final Pair<RecalDatumNode<T>, Double> nodeToRemove = getMinPenaltyNode();
        logger.info("Removing " + nodeToRemove.getFirst() + " with penalty " + nodeToRemove.getSecond());

        final Pair<RecalDatumNode<T>, Boolean> result = removeNode(nodeToRemove.getFirst());

        if ( ! result.getSecond() )
            throw new IllegalStateException("Never removed any node!");

        final RecalDatumNode<T> oneRemoved = result.getFirst();
        if ( oneRemoved == null )
            throw new IllegalStateException("Removed our root node, wow, didn't expect that");
        return oneRemoved;
    }

    private Pair<RecalDatumNode<T>, Double> getMinPenaltyNode() {
        final double myValue = isLeaf() ? Double.MAX_VALUE : getPenalty();
        Pair<RecalDatumNode<T>, Double> maxNode = new Pair<RecalDatumNode<T>, Double>(this, myValue);

        for ( final RecalDatumNode<T> sub : subnodes ) {
            final Pair<RecalDatumNode<T>, Double> subFind = sub.getMinPenaltyNode();
            if ( subFind.getSecond() < maxNode.getSecond() ) {
                maxNode = subFind;
            }
        }

        return maxNode;
    }

    private Pair<RecalDatumNode<T>, Boolean> removeNode(final RecalDatumNode<T> nodeToRemove) {
        if ( this == nodeToRemove ) {
            if ( isLeaf() )
                throw new IllegalStateException("Trying to remove a leaf node from the tree! " + this + " " + nodeToRemove);
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
                    final Pair<RecalDatumNode<T>, Boolean> maybeRemoved = sub1.removeNode(nodeToRemove);
                    removedSomething = maybeRemoved.getSecond();
                    sub.add(maybeRemoved.getFirst());
                }
            }

            final RecalDatumNode<T> node = new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, sub);
            return new Pair<RecalDatumNode<T>, Boolean>(node, removedSomething);
        }
    }
}
