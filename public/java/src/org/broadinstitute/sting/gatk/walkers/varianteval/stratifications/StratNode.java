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

package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Helper class representing a tree of stratification splits, where leaf nodes
 * are given a unique integer key starting at 0 and incrementing up to the
 * number of leaves in the tree.  This allows you to use this tree to produce
 * a key to map into an array index mapped data structure.
 *
 * Suppose I have to strats, each with two values: A = 1, 2 and B = 3, 4
 *
 * This data structure creates a tree such as:
 *
 * root -> A -> 1 -> B -> 3  : 0
 *                |- B -> 4  : 1
 *      |- A -> 2 -> B -> 3  : 2
 *                |- B -> 4  : 3
 *
 * This code allows us to efficiently look up a state key (A=2, B=3) and map it
 * to a specific key (an integer) that's unique over the tree
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
public class StratNode<T extends SetOfStates> implements Iterable<StratNode<T>> {
    int key = -1;
    final T stratifier;
    final Map<String, StratNode<T>> subnodes;

    public StratNode() {
        this.subnodes = Collections.emptyMap();
        this.stratifier = null;
    }

    StratNode(final T stratifier, final Map<String, StratNode<T>> subnodes) {
        this.stratifier = stratifier;
        this.subnodes = subnodes;
    }

    public void setKey(final int key) {
        if ( ! isLeaf() )
            throw new ReviewedStingException("Cannot set key of non-leaf node");
        this.key = key;
    }

    public int find(final List<String> states, int offset) {
        if ( isLeaf() ) // we're here!
            return key;
        else {
            final String state = states.get(offset);
            StratNode<T> subnode = subnodes.get(state);
            if ( subnode == null )
                throw new ReviewedStingException("Couldn't find state for " + state + " at node " + this);
            else
                return subnode.find(states, offset+1);
        }
    }

    public int getKey() {
        if ( ! isLeaf() )
            throw new ReviewedStingException("Cannot get key of non-leaf node");
        else
            return key;
    }

    protected Map<String, StratNode<T>> getSubnodes() {
        return subnodes;
    }

    public int size() {
        if ( isLeaf() )
            return 1;
        else {
            return subnodes.values().iterator().next().size() * subnodes.size();
        }
    }

    public T getSetOfStates() {
        return stratifier;
    }

    public boolean isLeaf() { return stratifier == null; }

    @Override
    public Iterator<StratNode<T>> iterator() {
        return new StratNodeIterator<T>(this);
    }
}
