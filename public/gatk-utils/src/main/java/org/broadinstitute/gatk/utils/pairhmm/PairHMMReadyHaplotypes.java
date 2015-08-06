/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.pairhmm;

import java.util.*;

/**
 * Collection of haplotypes sorted in a conveniently way to be run efficiently by the PairHMM.
 *
 * TODO not yet in use but likely to be as part of making graph-base likelihood run faster.
 * TODO this could be extended to the classical PairHMM implementation simplifyling the PairHMM API.
 */
public class PairHMMReadyHaplotypes implements Iterable<PairHMMReadyHaplotypes.Entry> {


    public class Entry {

        private final byte[] bases;

        private double likelihood = Double.NaN;

        protected Entry(final byte[] bases) {
            this.bases = bases;
        }

        public byte[] getBases() {
            return bases;
        }

        public void setLikelihood(final double lk) {
            likelihood = lk;
        }

        public double getLikelihood() {
            return likelihood;
        }

    }

    private Map<Entry,Map<Entry,Integer>> commonPrefixLength;

    private SortedSet<Entry> entries;

    private int capacity;

    private final Comparator<Entry> comparator = new Comparator<Entry>() {
        @Override
        public int compare(final Entry o1, final Entry o2) {
            final byte[] b1 = o1.bases;
            final byte[] b2 = o2.bases;
            Map<Entry,Integer> b1map = commonPrefixLength.get(o1);
            if (b1map == null)
                commonPrefixLength.put(o1, b1map = new HashMap<>(capacity));
            Map<Entry,Integer> b2map = commonPrefixLength.get(o2);
            if (b2map == null)
                commonPrefixLength.put(o2, b2map = new HashMap<>(capacity));
            final Integer previousI = b1map.get(o2) == null ? null : b1map.get(o2);
            int i;
            int result;
            final int iLimit = Math.min(b1.length,b2.length);
            if (previousI == null) {
                for (i = 0; i < iLimit; i++)
                    if (b1[i] != b2[i])
                        break;
                b1map.put(o2,i);
                b2map.put(o1,i);
            } else
                i = previousI;

            if (i < iLimit)
                result = Byte.compare(b1[i],b2[i]);
            else if (b1.length == b2.length)
                result = 0;
            else
                result = b1.length < b2.length ? -1 : 1;
            return result;
        }
    };

    public PairHMMReadyHaplotypes(final int capacity) {
        commonPrefixLength = new HashMap<>(capacity);
        entries = new TreeSet<>(comparator);
    }

    public void add(final byte[] bases) {
        final Entry entry = new Entry(bases);
        entries.add(entry);
    }

    public int size() {
        return entries.size();
    }

    @Override
    public Iterator iterator() {
        return new Iterator();
    }

    public class Iterator implements java.util.Iterator<Entry> {

        private java.util.Iterator<Entry> actualIterator;
        private Entry previousEntry;
        private Entry currentEntry;
        private int startIndex;
        private int cmp;

        private Iterator() {
            actualIterator = entries.iterator();
        }

        public boolean hasNext() {
            return actualIterator.hasNext();
        }

        public Entry next() {
            previousEntry = currentEntry;
            final Entry result = currentEntry = actualIterator.next();
            startIndex = -1;
            return result;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        public byte[] bases() {
            if (currentEntry == null)
                throw new NoSuchElementException();
            return currentEntry.bases;
        }

        public int startIndex() {
            if (startIndex >= 0)
                return startIndex;
            else if (previousEntry == null)
                return startIndex = 0;
            else {
                // The comparator will make sure the common-prefix-length is updated.
                // The result in a field so that we avoid dead code elimination.
                // perhaps I a bit paranohic but it does not harm to prevent.
                cmp = comparator.compare(previousEntry,currentEntry);
                return startIndex = commonPrefixLength.get(previousEntry).get(currentEntry);
            }
        }

        @Override
        public String toString() {
            return super.toString() + " cmp = " + cmp;
        }

        public void setLikelihood(final double likelihood) {
            if (currentEntry == null)
                throw new NoSuchElementException();
            currentEntry.setLikelihood(likelihood);
        }
    }

}
