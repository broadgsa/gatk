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

package org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

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
 * Note the structure of this tree is that the keys are -1 for all internal nodes, and
 * leafs are the only nodes with meaningful keys.  So for a tree with 2N nodes N of these
 * will be internal, with no keys, and meaningful maps from states -> subtrees.  The
 * other N nodes are leafs, with meaningful keys, empty maps, and null stratification objects
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
@Invariant({
        "(isLeaf() && stratifier == null && subnodes.isEmpty()) || (!isLeaf() && stratifier != null && !subnodes.isEmpty())"})
class StratNode<T extends Stratifier> implements Iterable<StratNode<T>> {
    int key = -1;
    final T stratifier;
    final Map<Object, StratNode<T>> subnodes; // NOTE, because we don't iterate our best option is a HashMap

    protected StratNode() {
        this.subnodes = Collections.emptyMap();
        this.stratifier = null;
    }

    protected StratNode(final T stratifier, final Map<Object, StratNode<T>> subnodes) {
        this.stratifier = stratifier;
        // important to reallocate an unmodififable hashmap with this specific size for space and safety
        this.subnodes = Collections.unmodifiableMap(new HashMap<Object, StratNode<T>>(subnodes));
    }

    @Requires("key >= 0")
    public void setKey(final int key) {
        if ( ! isLeaf() )
            throw new ReviewedGATKException("Cannot set key of non-leaf node");
        this.key = key;
    }

    @Requires({
            "states != null",
            "offset >= 0",
            "offset <= states.size()"
            })
    public int find(final List<Object> states, int offset) {
        if ( isLeaf() ) // we're here!
            return key;
        else {
            final Object state = states.get(offset);
            StratNode<T> subnode = subnodes.get(state);
            if ( subnode == null )
                return -1;
            else
                return subnode.find(states, offset+1);
        }
    }

    @Requires({
            "multipleStates != null",
            "offset >= 0",
            "offset <= multipleStates.size()",
            "keys != null",
            "offset == multipleStates.size() || multipleStates.get(offset) != null"})
    public void find(final List<List<Object>> multipleStates, final int offset, final HashSet<Integer> keys) {
        if ( isLeaf() ) // we're here!
            keys.add(key);
        else {
            for ( final Object state : multipleStates.get(offset) ) {
                // loop over all of the states at this offset
                final StratNode<T> subnode = subnodes.get(state);
                if ( subnode == null )
                    throw new ReviewedGATKException("Couldn't find state for " + state + " at node " + this);
                else
                    subnode.find(multipleStates, offset+1, keys);
            }
        }
    }

    @Ensures("result >= 0")
    public int getKey() {
        if ( ! isLeaf() )
            throw new ReviewedGATKException("Cannot get key of non-leaf node");
        else
            return key;
    }

    protected Map<Object, StratNode<T>> getSubnodes() {
        return subnodes;
    }

    @Ensures("result >= 0")
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

    /**
     * @return true if this node is a leaf
     */
    public boolean isLeaf() {
        return stratifier == null;
    }

    /**
     * Returns an iterator over this node and all subnodes including internal and leaf nodes
     * @return
     */
    @Override
    @Ensures("result != null")
    public Iterator<StratNode<T>> iterator() {
        return new StratNodeIterator<T>(this);
    }
}
