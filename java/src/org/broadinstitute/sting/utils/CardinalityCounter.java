package org.broadinstitute.sting.utils;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: fromer
 * Date: Aug 19, 2010
 * Time: 9:33:54 AM
 * To change this template use File | Settings | File Templates.
 */

/*
* CardinalityCounter object allows user to iterate over all assignment of arbitrary-cardinality variables.
 */
public class CardinalityCounter implements Iterator<int[]>, Iterable<int[]> {
    private int[] cards;
    private int[] valList;
    private boolean hasNext;

    public CardinalityCounter(int[] cards) {
        this.cards = cards;
        this.valList = new int[cards.length];
        for (int i = 0; i < cards.length; i++) {
            if (this.cards[i] <= 0)
                throw new IllegalArgumentException("CANNOT have zero cardinalities!");
            this.valList[i] = 0;
        }
        this.hasNext = true;
    }

    public boolean hasNext() {
        return hasNext;
    }

    public int[] next() {
        if (!hasNext())
            throw new StingException("CANNOT iterate past end!");

        // Copy the assignment to be returned:
        int[] nextList = new int[valList.length];
        for (int i = 0; i < valList.length; i++)
            nextList[i] = valList[i];

        // Find the assignment after this one:
        hasNext = false;
        int i = cards.length - 1;
        for (; i >= 0; i--) {
            if (valList[i] < (cards[i] - 1)) {
                valList[i]++;
                hasNext = true;
                break;
            }
            valList[i] = 0;
        }

        return nextList;
    }

    public void remove() {
        throw new RuntimeException("Cannot remove from CardinalityCounter!");
    }

    public Iterator<int[]> iterator() {
        return this;
    }
}
