package org.broadinstitute.sting.utils.recalibration;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

/**
 * A tree of recal datum, where each contains a set of sub datum representing sub-states of the higher level one
 *
 * @author Mark DePristo
 * @since 07/27/12
 */
public class RecalDatumNode<T extends RecalDatum> {
    private final static double SMALLEST_CHI2_PVALUE = 1e-300;
    protected static final Logger logger = Logger.getLogger(RecalDatumNode.class);

    /**
     * fixedPenalty is this value if it's considered fixed
     */
    private final static double UNINITIALIZED = Double.NEGATIVE_INFINITY;

    private final T recalDatum;
    private double fixedPenalty = UNINITIALIZED;
    private final Set<RecalDatumNode<T>> subnodes;

    @Requires({"recalDatum != null"})
    public RecalDatumNode(final T recalDatum) {
        this(recalDatum, new HashSet<RecalDatumNode<T>>());
    }

    @Override
    public String toString() {
        return recalDatum.toString();
    }

    @Requires({"recalDatum != null", "subnodes != null"})
    public RecalDatumNode(final T recalDatum, final Set<RecalDatumNode<T>> subnodes) {
        this(recalDatum, UNINITIALIZED, subnodes);
    }

    @Requires({"recalDatum != null"})
    protected RecalDatumNode(final T recalDatum, final double fixedPenalty) {
        this(recalDatum, fixedPenalty, new HashSet<RecalDatumNode<T>>());
    }

    @Requires({"recalDatum != null", "subnodes != null"})
    protected RecalDatumNode(final T recalDatum, final double fixedPenalty, final Set<RecalDatumNode<T>> subnodes) {
        this.recalDatum = recalDatum;
        this.fixedPenalty = fixedPenalty;
        this.subnodes = new HashSet<RecalDatumNode<T>>(subnodes);
    }

    /**
     * Get the recal data associated with this node
     * @return
     */
    @Ensures("result != null")
    public T getRecalDatum() {
        return recalDatum;
    }

    /**
     * The set of all subnodes of this tree.  May be modified.
     * @return
     */
    @Ensures("result != null")
    public Set<RecalDatumNode<T>> getSubnodes() {
        return subnodes;
    }

    /**
     * Return the fixed penalty, if set, or else the the calculated penalty for this node
     * @return
     */
    public double getPenalty() {
        if ( fixedPenalty != UNINITIALIZED )
            return fixedPenalty;
        else
            return calcPenalty();
    }

    /**
     * Set the fixed penalty for this node to a fresh calculation from calcPenalty
     *
     * This is important in the case where you want to compute the penalty from a full
     * tree and then chop the tree up afterwards while considering the previous penalties.
     * If you don't call this function then manipulating the tree may result in the
     * penalty functions changing with changes in the tree.
     *
     * @param doEntireTree recurse into all subnodes?
     * @return the fixed penalty for this node
     */
    public double calcAndSetFixedPenalty(final boolean doEntireTree) {
        fixedPenalty = calcPenalty();
        if ( doEntireTree )
            for ( final RecalDatumNode<T> sub : subnodes )
                sub.calcAndSetFixedPenalty(doEntireTree);
        return fixedPenalty;
    }

    /**
     * Add node to the set of subnodes of this node
     * @param sub
     */
    @Requires("sub != null")
    public void addSubnode(final RecalDatumNode<T> sub) {
        subnodes.add(sub);
    }

    /**
     * Is this a leaf node (i.e., has no subnodes)?
     * @return
     */
    public boolean isLeaf() {
        return subnodes.isEmpty();
    }

    /**
     * Is this node immediately above only leaf nodes?
     *
     * @return
     */
    public boolean isAboveOnlyLeaves() {
        for ( final RecalDatumNode<T> sub : subnodes )
            if ( ! sub.isLeaf() )
                return false;
        return true;
    }

    /**
     * What's the immediate number of subnodes from this node?
     * @return
     */
    @Ensures("result >= 0")
    public int getNumSubnodes() {
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

    /**
     * The maximum penalty among all nodes
     * @return
     */
    public double maxPenalty(final boolean leafOnly) {
        double max = ! leafOnly || isLeaf() ? getPenalty() : Double.MIN_VALUE;
        for ( final RecalDatumNode<T> sub : subnodes )
            max = Math.max(max, sub.maxPenalty(leafOnly));
        return max;
    }

    /**
     * The minimum penalty among all nodes
     * @return
     */
    public double minPenalty(final boolean leafOnly) {
        double min = ! leafOnly || isLeaf() ? getPenalty() : Double.MAX_VALUE;
        for ( final RecalDatumNode<T> sub : subnodes )
            min = Math.min(min, sub.minPenalty(leafOnly));
        return min;
    }

    /**
     * What's the longest branch from this node to any leaf?
     * @return
     */
    public int maxDepth() {
        int subMax = 0;
        for ( final RecalDatumNode<T> sub : subnodes )
            subMax = Math.max(subMax, sub.maxDepth());
        return subMax + 1;
    }

    /**
     * What's the shortest branch from this node to any leaf?  Includes this node
     * @return
     */
    @Ensures("result > 0")
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

    /**
     * Return the number of nodes, including this one, reachable from this node
     * @return
     */
    @Ensures("result > 0")
    public int size() {
        int size = 1;
        for ( final RecalDatumNode<T> sub : subnodes )
            size += sub.size();
        return size;
    }

    /**
     * Count the number of leaf nodes reachable from this node
     *
     * @return
     */
    @Ensures("result >= 0")
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

    /**
     * Calculate the phred-scaled p-value for a chi^2 test for independent among subnodes of this node.
     *
     * The chi^2 value indicates the degree of independence of the implied error rates among the
     * immediate subnodes
     *
     * @return the phred-scaled p-value for chi2 penalty, or 0.0 if it cannot be calculated
     */
    private double calcPenalty() {
        if ( isLeaf() || freeToMerge() )
            return 0.0;
        else if ( subnodes.size() == 1 )
            // only one value, so its free to merge away
            return 0.0;
        else {
            final long[][] counts = new long[subnodes.size()][2];

            int i = 0;
            for ( final RecalDatumNode<T> subnode : subnodes ) {
                // use the yates correction to help avoid all zeros => NaN
                counts[i][0] = Math.round(subnode.getRecalDatum().getNumMismatches()) + 1L;
                counts[i][1] = Math.round(subnode.getRecalDatum().getNumObservations()) + 2L;
                i++;
            }

            try {
                final double chi2PValue = new ChiSquareTestImpl().chiSquareTest(counts);
                final double penalty = -10.0 * Math.log10(Math.max(chi2PValue, SMALLEST_CHI2_PVALUE));

                // make sure things are reasonable and fail early if not
                if (Double.isInfinite(penalty) || Double.isNaN(penalty))
                    throw new ReviewedStingException("chi2 value is " + chi2PValue + " at " + getRecalDatum());

                return penalty;
            } catch ( MathException e ) {
                throw new ReviewedStingException("Failed in calculating chi2 value", e);
            }
        }
    }

    /**
     * Is this node free to merge because its rounded Q score is the same as all nodes below
     * @return
     */
    private boolean freeToMerge() {
        if ( isLeaf() ) // leaves are free to merge
            return true;
        else {
            final byte myQual = getRecalDatum().getEmpiricalQualityAsByte();
            for ( final RecalDatumNode<T> sub : subnodes )
                if ( sub.getRecalDatum().getEmpiricalQualityAsByte() != myQual )
                    return false;
            return true;
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

    /**
     * Return a freshly allocated tree prunes to have no more than maxDepth from the root to any leaf
     *
     * @param maxDepth
     * @return
     */
    public RecalDatumNode<T> pruneToDepth(final int maxDepth) {
        if ( maxDepth < 1 )
            throw new IllegalArgumentException("maxDepth < 1");
        else {
            final Set<RecalDatumNode<T>> subPruned = new HashSet<RecalDatumNode<T>>(getNumSubnodes());
            if ( maxDepth > 1 )
                for ( final RecalDatumNode<T> sub : subnodes )
                    subPruned.add(sub.pruneToDepth(maxDepth - 1));
            return new RecalDatumNode<T>(getRecalDatum(), fixedPenalty, subPruned);
        }
    }

    /**
     * Return a freshly allocated tree with to no more than maxElements in order of penalty
     *
     * Note that nodes must have fixed penalties to this algorithm will fail.
     *
     * @param maxElements
     * @return
     */
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
     * Return a freshly allocated tree where all mergable nodes with < maxPenalty are merged
     *
     * Note that nodes must have fixed penalties to this algorithm will fail.
     *
     * @param maxPenaltyIn the maximum penalty we are allowed to incur for a merge
     * @param applyBonferroniCorrection if true, we will adjust penalty by the phred-scaled bonferroni correction
     *                                  for the size of the initial tree.  That is, if there are 10 nodes in the
     *                                  tree and maxPenalty is 20 we will actually enforce 10^-2 / 10 = 10^-3 = 30
     *                                  penalty for multiple testing
     * @return
     */
    public RecalDatumNode<T> pruneToNoMoreThanPenalty(final double maxPenaltyIn, final boolean applyBonferroniCorrection) {
        RecalDatumNode<T> root = this;

        final double bonferroniCorrection = 10 * Math.log10(this.size());
        final double maxPenalty = applyBonferroniCorrection ? maxPenaltyIn + bonferroniCorrection : maxPenaltyIn;

        if ( applyBonferroniCorrection )
        logger.info(String.format("Applying Bonferroni correction for %d nodes = %.2f to initial penalty %.2f for total " +
                "corrected max penalty of %.2f", this.size(), bonferroniCorrection, maxPenaltyIn, maxPenalty));

        while ( true ) {
            final Pair<RecalDatumNode<T>, Double> minPenaltyNode = root.getMinPenaltyAboveLeafNode();

            if ( minPenaltyNode == null || minPenaltyNode.getSecond() > maxPenalty ) {
                // nothing to merge, or the best candidate is above our max allowed
                if ( minPenaltyNode == null ) {
                    if ( logger.isDebugEnabled() ) logger.debug("Stopping because no candidates could be found");
                } else {
                    if ( logger.isDebugEnabled() ) logger.debug("Stopping because node " + minPenaltyNode.getFirst() + " has penalty " + minPenaltyNode.getSecond() + " > max " + maxPenalty);
                }
                break;
            } else {
                // remove the lowest penalty element, and continue
                if ( logger.isDebugEnabled() ) logger.debug("Removing node " + minPenaltyNode.getFirst() + " with penalty " + minPenaltyNode.getSecond());
                root = root.removeLowestPenaltyNode();
            }
        }

        // no more candidates exist with penalty < maxPenalty
        return root;
    }


    /**
     * Find the lowest penalty above leaf node in the tree, and return a tree without it
     *
     * Note this excludes the current (root) node
     *
     * @return
     */
    private RecalDatumNode<T> removeLowestPenaltyNode() {
        final Pair<RecalDatumNode<T>, Double> nodeToRemove = getMinPenaltyAboveLeafNode();
        if ( logger.isDebugEnabled() )
            logger.debug("Removing " + nodeToRemove.getFirst() + " with penalty " + nodeToRemove.getSecond());

        final Pair<RecalDatumNode<T>, Boolean> result = removeNode(nodeToRemove.getFirst());

        if ( ! result.getSecond() )
            throw new IllegalStateException("Never removed any node!");

        final RecalDatumNode<T> oneRemoved = result.getFirst();
        if ( oneRemoved == null )
            throw new IllegalStateException("Removed our root node, wow, didn't expect that");
        return oneRemoved;
    }

    /**
     * Finds in the tree the node with the lowest penalty whose subnodes are all leaves
     *
     * @return the node and its penalty, or null if no such node exists
     */
    private Pair<RecalDatumNode<T>, Double> getMinPenaltyAboveLeafNode() {
        if ( isLeaf() )
            // not allowed to remove leafs directly
            return null;
        if ( isAboveOnlyLeaves() )
            // we only consider removing nodes above all leaves
            return new Pair<RecalDatumNode<T>, Double>(this, getPenalty());
        else {
            // just recurse, taking the result with the min penalty of all subnodes
            Pair<RecalDatumNode<T>, Double> minNode = null;
            for ( final RecalDatumNode<T> sub : subnodes ) {
                final Pair<RecalDatumNode<T>, Double> subFind = sub.getMinPenaltyAboveLeafNode();
                if ( subFind != null && (minNode == null || subFind.getSecond() < minNode.getSecond()) ) {
                    minNode = subFind;
                }
            }
            return minNode;
        }
    }

    /**
     * Return a freshly allocated tree without the node nodeToRemove
     *
     * @param nodeToRemove
     * @return
     */
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
            final Set<RecalDatumNode<T>> sub = new HashSet<RecalDatumNode<T>>(getNumSubnodes());

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

    /**
     * Return a collection of all of the data in the leaf nodes of this tree
     *
     * @return
     */
    public Collection<T> getAllLeaves() {
        final LinkedList<T> list = new LinkedList<T>();
        getAllLeavesRec(list);
        return list;
    }

    /**
     * Helpful recursive function for getAllLeaves()
     *
     * @param list the destination for the list of leaves
     */
    private void getAllLeavesRec(final LinkedList<T> list) {
        if ( isLeaf() )
            list.add(getRecalDatum());
        else {
            for ( final RecalDatumNode<T> sub : subnodes )
                sub.getAllLeavesRec(list);
        }
    }
}
