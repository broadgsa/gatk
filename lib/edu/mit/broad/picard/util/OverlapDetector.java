package edu.mit.broad.picard.util;

import java.util.*;

/**
 * Utility class to efficiently do in memory overlap detection between a large
 * set of mapping like objects, and one or more candidate mappings.
 */
public class OverlapDetector<T> {
    private Map<Object, IntervalTree<Set<T>>> cache = new HashMap<Object, IntervalTree<Set<T>>>();
    private final int lhsBuffer;
    private final int rhsBuffer;

    /**
     * Constructs an overlap detector.
     * @param lhsBuffer the amount by which to "trim" coordinates of mappings on the left
     *                  hand side when calculating overlaps
     * @param rhsBuffer the amount by which to "trim" coordinates of mappings on the right
     *                  hand side when calculating overlaps
     */
    public OverlapDetector(int lhsBuffer, int rhsBuffer) {
        this.lhsBuffer = lhsBuffer;
        this.rhsBuffer = rhsBuffer;
    }

    /** Adds a mapping to the set of mappings against which to match candidates. */
    public void addLhs(T object, Interval interval) {
        Object seqId = interval.getSequence();

        IntervalTree<Set<T>> tree = this.cache.get(seqId);
        if (tree == null) {
            tree = new IntervalTree<Set<T>>();
            this.cache.put(seqId, tree);
        }

        int start = interval.getStart() + this.lhsBuffer;
        int end   = interval.getEnd()   - this.lhsBuffer;

        Set<T> objects = new HashSet<T>();
        objects.add(object);
        if (start <= end)  // Don't put in sequences that have no overlappable bases
        {
            Set<T> alreadyThere = tree.put(start, end, objects);
            if (alreadyThere != null)
            {
                alreadyThere.add(object);
                tree.put(start, end, alreadyThere);
            }
        }
    }

    /** Adds all items to the overlap detector. */
    public void addAll(List<T> objects, List<Interval> intervals) {
        if (objects.size() != intervals.size()) {
            throw new IllegalArgumentException("Objects and intervals must be the same size.");
        }

        for (int i=0; i<objects.size(); ++i) {
            addLhs(objects.get(i), intervals.get(i));
        }
    }

    /** Gets the collection of objects that overlap the provided mapping. */
    public Collection<T> getOverlaps(Interval rhs)  {
        Collection<T> matches = new ArrayList<T>();

        Object seqId = rhs.getSequence();
        IntervalTree<Set<T>> tree = this.cache.get(seqId);
        int start = rhs.getStart() + this.rhsBuffer;
        int end = rhs.getEnd() - this.rhsBuffer;

        if (tree != null && start <= end)
        {
            Iterator<IntervalTree.Node<Set<T>>> it = tree.overlappers(start, end);
            while (it.hasNext())
            {
                IntervalTree.Node<Set<T>> node = it.next();
                matches.addAll(node.getValue());
            }
        }

        return matches;
    }

    /** Gets all the objects that could be returned by the overlap detector. */
    public Collection<T> getAll() {
        Collection<T> all = new HashSet<T>();
        for (IntervalTree<Set<T>> tree : this.cache.values()) {
            for (IntervalTree.Node<Set<T>> node : tree) {
                all.addAll(node.getValue());
            }
        }

        return all;
    }
}
