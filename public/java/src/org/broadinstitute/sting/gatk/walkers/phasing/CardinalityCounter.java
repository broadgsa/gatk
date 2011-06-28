/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Iterator;

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
            throw new ReviewedStingException("CANNOT iterate past end!");

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
