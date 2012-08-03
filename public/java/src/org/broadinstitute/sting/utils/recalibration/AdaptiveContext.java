package org.broadinstitute.sting.utils.recalibration;

import java.util.*;

/**
 * Functions for working with AdaptiveContexts
 *
 * User: depristo
 * Date: 8/3/12
 * Time: 12:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class AdaptiveContext {
    private AdaptiveContext() {}

    /**
     * Return a freshly allocated tree filled in completely to fillDepth with
     * all combinations of {A,C,G,T}^filldepth contexts.  For nodes
     * in the tree, they are simply copied.  When the algorithm needs to
     * generate new nodes (because they are missing) the subnodes inherit the
     * observation and error counts of their parent.
     *
     * This algorithm produces data consistent with the standard output in a BQSR recal
     * file for the Context covariate
     *
     * @param root
     * @param fillDepth
     * @return
     */
    public static RecalDatumNode<ContextDatum> fillToDepth(final RecalDatumNode<ContextDatum> root, final int fillDepth) {
        if ( root == null ) throw new IllegalArgumentException("root is null");
        if ( fillDepth < 0 ) throw new IllegalArgumentException("fillDepth is < 0");

        return fillToDepthRec(root, fillDepth, 0);
    }

    private static RecalDatumNode<ContextDatum> fillToDepthRec(final RecalDatumNode<ContextDatum> parent,
                                                              final int fillDepth,
                                                              final int currentDepth) {
        // three cases:
        //   We are in the tree and so just recursively build
        //   We have reached our depth goal, so just return the parent since we are done
        //   We are outside of the tree, in which case we need to pointer to our parent node so we can
        //     we info (N, M) and we need a running context
        if ( currentDepth < fillDepth ) {
            // we need to create subnodes for each base, and propogate N and M down
            final RecalDatumNode<ContextDatum> newParent = new RecalDatumNode<ContextDatum>(parent.getRecalDatum());

            for ( final String base : Arrays.asList("A", "C", "G", "T")) {
                ContextDatum subContext;
                Set<RecalDatumNode<ContextDatum>> subContexts;

                final RecalDatumNode<ContextDatum> subNode = findSubcontext(parent.getRecalDatum().context + base, parent);
                if ( subNode != null ) {
                    // we have a subnode corresponding to the expected one, just copy and recurse
                    subContext = subNode.getRecalDatum();
                    subContexts = subNode.getSubnodes();
                } else {
                    // have to create a new one
                    subContext = new ContextDatum(parent.getRecalDatum().context + base,
                            parent.getRecalDatum().getNumObservations(), parent.getRecalDatum().getNumMismatches());
                    subContexts = Collections.emptySet();
                }

                newParent.addSubnode(
                        fillToDepthRec(new RecalDatumNode<ContextDatum>(subContext, subContexts),
                                fillDepth, currentDepth + 1));
            }
            return newParent;
        } else {
            return parent;
        }
    }

    /**
     * Go from a flat list of contexts to the tree implied by the contexts
     *
     * Implicit nodes are created as needed, and their observation and error counts are the sum of the
     * all of their subnodes.
     *
     * Note this does not guarentee the tree is complete, as some contexts (e.g., AAT) may be missing
     * from the tree because they are absent from the input list of contexts.
     *
     * For input AAG, AAT, AC, G would produce the following tree:
     *
     * - x [root]
     *   - A
     *     - A
     *       - T
     *       - G
     *     - C
     *   - G
     *
     * sets the fixed penalties in the resulting tree as well
     *
     * @param flatContexts list of flat contexts
     * @return
     */
    public static RecalDatumNode<ContextDatum> createTreeFromFlatContexts(final List<ContextDatum> flatContexts) {
        if ( flatContexts == null || flatContexts.isEmpty() )
            throw new IllegalArgumentException("flatContexts cannot be empty or null");

        final Queue<RecalDatumNode<ContextDatum>> remaining = new LinkedList<RecalDatumNode<ContextDatum>>();
        final Map<String, RecalDatumNode<ContextDatum>> contextToNodes = new HashMap<String, RecalDatumNode<ContextDatum>>();
        RecalDatumNode<ContextDatum> root = null;

        // initialize -- start with all of the contexts
        for ( final ContextDatum cd : flatContexts )
            remaining.add(new RecalDatumNode<ContextDatum>(cd));

        while ( remaining.peek() != null ) {
            final RecalDatumNode<ContextDatum> add = remaining.poll();
            final ContextDatum cd = add.getRecalDatum();

            final String parentContext = cd.getParentContext();
            RecalDatumNode<ContextDatum> parent = contextToNodes.get(parentContext);
            if ( parent == null ) {
                // haven't yet found parent, so make one, and enqueue it for processing
                parent = new RecalDatumNode<ContextDatum>(new ContextDatum(parentContext, 0, 0));
                contextToNodes.put(parentContext, parent);

                if ( parentContext != ContextDatum.ROOT_CONTEXT )
                    remaining.add(parent);
                else
                    root = parent;
            }

            parent.getRecalDatum().incrementNumObservations(cd.getNumObservations());
            parent.getRecalDatum().incrementNumMismatches(cd.getNumMismatches());
            parent.addSubnode(add);
        }

        if ( root == null )
            throw new RuntimeException("root is unexpectedly null");

        // set the fixed penalty everywhere in the tree, so that future modifications don't change the penalties
        root.calcAndSetFixedPenalty(true);

        return root;
    }

    /**
     * Finds immediate subnode with contextToFind, or null if none exists
     *
     * @param tree whose subnodes should be searched
     * @return
     */
    public static RecalDatumNode<ContextDatum> findSubcontext(final String contextToFind, final RecalDatumNode<ContextDatum> tree) {
        for ( final RecalDatumNode<ContextDatum> sub : tree.getSubnodes() )
            if ( sub.getRecalDatum().context.equals(contextToFind) )
                return sub;
        return null;
    }
}
