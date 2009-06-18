package org.broadinstitute.sting.playground.piecemealannotator;

import java.util.PriorityQueue;

public class BoundedScoringSet<E extends Comparable<E> > {
    private PriorityQueue<E> pq;
    private int maximumSize;

    public BoundedScoringSet(int maximumSize) {
        pq = new PriorityQueue<E>(maximumSize);
        this.maximumSize = maximumSize;
    }

    public boolean add(E o) {
        if (canAdd(o)) {
            pq.add(o);

            while (pq.size() > maximumSize) {
                pq.poll();
            }

            return true;
        }

        return false;
    }

    private boolean canAdd(E o) { return pq.size() < maximumSize || o.compareTo(pq.peek()) == 1; }

    public void clear() { pq.clear(); }

    public boolean contains(E o) { return pq.contains(o); }

    public boolean offer(E o) { return pq.offer(o); }

    public E peek() { return pq.peek(); }

    public E poll() { return pq.poll(); }

    public boolean remove(E o) { return pq.remove(o); }

    public int size() { return pq.size(); }

    public E[] toArray(E[] os) { return pq.toArray(os); }
}
