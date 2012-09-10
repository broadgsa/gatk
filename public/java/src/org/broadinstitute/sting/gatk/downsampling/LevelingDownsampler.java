/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.downsampling;

import org.broadinstitute.sting.utils.MathUtils;

import java.util.*;

/**
 * Leveling Downsampler: Given a set of Lists of arbitrary items and a target size, removes items from
 * the Lists in an even fashion until the total size of all Lists is <= the target size. Leveling
 * does not occur until all Lists have been submitted and signalEndOfInput() is called.
 *
 * The Lists should be LinkedLists for maximum efficiency during item removal, however other
 * kinds of Lists are also accepted (albeit at a slight performance penalty).
 *
 * Since this downsampler extends the Downsampler interface rather than the ReadsDownsampler interface,
 * the Lists need not contain reads. However this downsampler may not be wrapped within one of the
 * DownsamplingReadsIterators
 *
 * @param <T> the List type representing the stacks to be leveled
 * @param <E> the type of the elements of each List
 *
 * @author David Roazen
 */
public class LevelingDownsampler<T extends List<E>, E> implements Downsampler<T> {

    private int targetSize;

    private List<T> groups;

    private boolean groupsAreFinalized;

    private int numDiscardedItems;

    /**
     * Construct a LevelingDownsampler
     *
     * @param targetSize the sum of the sizes of all individual Lists this downsampler is fed may not exceed
     *                   this value -- if it does, items are removed from Lists evenly until the total size
     *                   is <= this value
     */
    public LevelingDownsampler( int targetSize ) {
        this.targetSize = targetSize;
        clear();
        reset();
    }

    public void submit( T item ) {
        groups.add(item);
    }

    public void submit( Collection<T> items ){
        groups.addAll(items);
    }

    public boolean hasFinalizedItems() {
        return groupsAreFinalized && groups.size() > 0;
    }

    public List<T> consumeFinalizedItems() {
        if ( ! hasFinalizedItems() ) {
            return new ArrayList<T>();
        }

        // pass by reference rather than make a copy, for speed
        List<T> toReturn = groups;
        clear();
        return toReturn;
    }

    public boolean hasPendingItems() {
        return ! groupsAreFinalized && groups.size() > 0;
    }

    public T peekFinalized() {
        return hasFinalizedItems() ? groups.get(0) : null;
    }

    public T peekPending() {
        return hasPendingItems() ? groups.get(0) : null;
    }

    public int getNumberOfDiscardedItems() {
        return numDiscardedItems;
    }

    public void signalEndOfInput() {
        levelGroups();
        groupsAreFinalized = true;
    }

    public void clear() {
        groups = new ArrayList<T>();
        groupsAreFinalized = false;
    }

    public void reset() {
        numDiscardedItems = 0;
    }

    private void levelGroups() {
        int totalSize = 0;
        int[] groupSizes = new int[groups.size()];
        int currentGroupIndex = 0;

        for ( T group : groups ) {
            groupSizes[currentGroupIndex] = group.size();
            totalSize += groupSizes[currentGroupIndex];
            currentGroupIndex++;
        }

        if ( totalSize <= targetSize ) {
            return;    // no need to eliminate any items
        }

        // We will try to remove exactly this many items, however we will refuse to allow any
        // one group to fall below size 1, and so might end up removing fewer items than this
        int numItemsToRemove = totalSize - targetSize;

        currentGroupIndex = 0;
        int numConsecutiveUmodifiableGroups = 0;

        // Continue until we've either removed all the items we wanted to, or we can't
        // remove any more items without violating the constraint that all groups must
        // be left with at least one item
        while ( numItemsToRemove > 0 && numConsecutiveUmodifiableGroups < groupSizes.length ) {
            if ( groupSizes[currentGroupIndex] > 1 ) {
                groupSizes[currentGroupIndex]--;
                numItemsToRemove--;
                numConsecutiveUmodifiableGroups = 0;
            }
            else {
                numConsecutiveUmodifiableGroups++;
            }

            currentGroupIndex = (currentGroupIndex + 1) % groupSizes.length;
        }

        // Now we actually go through and reduce each group to its new count as specified in groupSizes
        currentGroupIndex = 0;
        for ( T group : groups ) {
            downsampleOneGroup(group, groupSizes[currentGroupIndex]);
            currentGroupIndex++;
        }
    }

    private void downsampleOneGroup( T group, int numItemsToKeep ) {
        if ( numItemsToKeep >= group.size() ) {
            return;
        }

        numDiscardedItems += group.size() - numItemsToKeep;

        BitSet itemsToKeep = new BitSet(group.size());
        for ( Integer selectedIndex : MathUtils.sampleIndicesWithoutReplacement(group.size(), numItemsToKeep) ) {
            itemsToKeep.set(selectedIndex);
        }

        int currentIndex = 0;

        // If our group is a linked list, we can remove the desired items in a single O(n) pass with an iterator
        if ( group instanceof LinkedList ) {
            Iterator iter = group.iterator();
            while ( iter.hasNext() ) {
                iter.next();

                if ( ! itemsToKeep.get(currentIndex) ) {
                    iter.remove();
                }

                currentIndex++;
            }
        }
        // If it's not a linked list, it's more efficient to copy the desired items into a new list and back rather
        // than suffer O(n^2) of item shifting
        else {
            List<E> keptItems = new ArrayList<E>(numItemsToKeep);

            for ( E item : group ) {
                if ( itemsToKeep.get(currentIndex) ) {
                    keptItems.add(item);
                }
                currentIndex++;
            }
            group.clear();
            group.addAll(keptItems);
        }
    }
}
