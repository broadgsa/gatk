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

import java.util.*;

/**
 * Represents the full state space of all stratification combinations
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
public class StratificationStates<T extends SetOfStates> {
    private final StratNode<T> root;

    public StratificationStates(final List<T> strats) {
        this.root = buildStratificationTree(new LinkedList<T>(strats));

        assignKeys(root, 0);
    }

    private StratNode<T> buildStratificationTree(final Queue<T> strats) {
        final T first = strats.poll();
        if ( first == null ) {
            // we are at a leaf
            return new StratNode<T>();
        } else {
            // we are in the middle of the tree
            final Collection<Object> states = first.getAllStates();
            final LinkedHashMap<Object, StratNode<T>> subNodes = new LinkedHashMap<Object, StratNode<T>>(states.size());
            for ( final Object state : states ) {
                // have to copy because poll modifies the queue
                final Queue<T> copy = new LinkedList<T>(strats);
                subNodes.put(state, buildStratificationTree(copy));
            }
            return new StratNode<T>(first, subNodes);
        }
    }

    public int getNStates() {
        return root.size();
    }

    public StratNode<T> getRoot() {
        return root;
    }
    
    public int getKey(final List<Object> states) {
        return root.find(states, 0);
    }

    public Set<Integer> getKeys(final List<List<Object>> allStates) {
        final HashSet<Integer> keys = new HashSet<Integer>();
        root.find(allStates, 0, keys);
        return keys;
    }
    
    private void assignKeys(final StratNode<T> root, int key) {
        for ( final StratNode<T> node : root ) {
            if ( node.isLeaf() )
                node.setKey(key++);
        }
    }
    
    public static List<List<Object>> combineStates(final List<Object> first, final List<Object> second) {
        List<List<Object>> combined = new ArrayList<List<Object>>(first.size());
        for ( int i = 0; i < first.size(); i++ ) {
            final Object firstI = first.get(i);
            final Object secondI = second.get(i);
            if ( firstI.equals(secondI) ) 
                combined.add(Collections.singletonList(firstI));
            else 
                combined.add(Arrays.asList(firstI, secondI));
        }
        return combined;
    }
}
