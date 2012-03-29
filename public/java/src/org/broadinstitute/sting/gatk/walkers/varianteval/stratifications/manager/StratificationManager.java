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

package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Represents the full state space of all stratification combinations
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
public class StratificationManager<K extends SetOfStates, V> implements Map<List<Object>, V> {
    private final StratNode<K> root;
    private final int size;

    private final ArrayList<K> stratifiers;

    // values associated with each key
    private final ArrayList<V> valuesByKey;
    private final ArrayList<List<Object>> stratifierValuesByKey;

    // -------------------------------------------------------------------------------------
    //
    // creating the manager
    //
    // -------------------------------------------------------------------------------------

    @Requires("!strats.isEmpty()")
    public StratificationManager(final List<K> strats) {
        stratifiers = new ArrayList<K>(strats);
        this.root = buildStratificationTree(new LinkedList<K>(strats));
        assignKeys(root);

        this.size = root.size();
        if ( this.size == 0 )
            throw new ReviewedStingException("Size == 0 in StratificationManager");

        this.valuesByKey = new ArrayList<V>(size());
        this.stratifierValuesByKey = new ArrayList<List<Object>>(size());
        for ( int i = 0; i < size(); i++ ) {
            this.valuesByKey.add(null);
            this.stratifierValuesByKey.add(null);
        }
        assignStratifierValuesByKey(root);
    }

    private StratNode<K> buildStratificationTree(final Queue<K> strats) {
        final K first = strats.poll();
        if ( first == null ) {
            // we are at a leaf
            return new StratNode<K>();
        } else {
            // we are in the middle of the tree
            final Collection<Object> states = first.getAllStates();
            
            if ( states.isEmpty() )
                throw new ReviewedStingException("State " + first + " is empty!");
            
            final LinkedHashMap<Object, StratNode<K>> subNodes = new LinkedHashMap<Object, StratNode<K>>(states.size());
            for ( final Object state : states ) {
                // have to copy because poll modifies the queue
                final Queue<K> copy = new LinkedList<K>(strats);
                subNodes.put(state, buildStratificationTree(copy));
            }
            return new StratNode<K>(first, subNodes);
        }
    }

    @Requires("root == this.root")
    private void assignKeys(final StratNode<K> root) {
        int key = 0;
        for ( final StratNode<K> node : root ) {
            if ( node.isLeaf() )
                node.setKey(key++);
        }
    }

    public void assignStratifierValuesByKey(final StratNode<K> root) {
        assignStratifierValuesByKey(root, new LinkedList<Object>());
        
        for ( List<Object> stateValues : stratifierValuesByKey )
            if ( stateValues == null )
                throw new ReviewedStingException("Found a null state value set that's null");
    }

    public void assignStratifierValuesByKey(final StratNode<K> node, final LinkedList<Object> states) {
        if ( node.isLeaf() ) { // we're here!
            if ( states.isEmpty() )
                throw new ReviewedStingException("Found a leaf node with an empty state values vector");
            stratifierValuesByKey.set(node.getKey(), new ArrayList<Object>(states));
        } else {
            for ( Map.Entry<Object, StratNode<K>> entry : node.getSubnodes().entrySet() ) {
                final LinkedList<Object> newStates = new LinkedList<Object>(states);
                newStates.addLast(entry.getKey());
                assignStratifierValuesByKey(entry.getValue(), newStates);
            }
        }
    }
    
    // -------------------------------------------------------------------------------------
    //
    // simple accessors
    //
    // -------------------------------------------------------------------------------------

    @Ensures("result >= 0")
    public int size() {
        return size;
    }

    @Ensures("result != null")
    public StratNode<K> getRoot() {
        return root;
    }

    // -------------------------------------------------------------------------------------
    //
    // mapping from states -> keys
    //
    // -------------------------------------------------------------------------------------

    @Requires("states != null")
    @Ensures("result >= -1")
    public int getKey(final List<Object> states) {
        return root.find(states, 0);
    }

    @Requires("allStates != null")
    @Ensures("result != null")
    public Set<Integer> getKeys(final List<List<Object>> allStates) {
        final HashSet<Integer> keys = new HashSet<Integer>();
        root.find(allStates, 0, keys);
        return keys;
    }

    public Map<K, Object> getStateForKey(final int key) {
        final Map<K, Object> states = new HashMap<K, Object>(stratifiers.size());
        for ( int i = 0; i < stratifiers.size(); i++ ) {
            final K strat = stratifiers.get(i);
            final Object stratValue = stratifierValuesByKey.get(key).get(i);
            states.put(strat, stratValue);
        }
        return states;
    }

    // -------------------------------------------------------------------------------------
    //
    // valuesByKey
    //
    // -------------------------------------------------------------------------------------

    @Override
    @Ensures("result != null")
    public ArrayList<V> values() {
        return valuesByKey;
    }
    
    public Collection<V> values(List<List<Object>> states) {
        // TODO -- SHOULD BE INLINE TO AVOID CREATING LIST OF KEYS JUST TO ITERATE OVER IT
        Collection<V> vals = new LinkedList<V>();
        for ( int key : getKeys(states) ) 
            vals.add(get(key));
        return vals;
    }

    @Requires("key >= 0 && key <= size()")
    @Ensures("get(key) == value")
    public void set(final int key, final V value) {
        valuesByKey.set(key, value);
    }

    @Requires("key >= 0 && key <= size()")
    public V get(final int key) {
        return valuesByKey.get(key);
    }

    @Requires("getKey(states) != -1")
    public V get(final List<Object> states) {
        return get(getKey(states));
    }

    @Override
    public V get(final Object o) {
        return get((List<Object>)o);
    }

    @Override
    public boolean isEmpty() {
        return false;
    }

    public boolean containsKey(final List<Object> o) {
        return getKey(o) != -1;
    }

    @Override
    public boolean containsKey(final Object o) {
        return containsKey((List<Object>)o);
    }

    @Override
    public boolean containsValue(final Object o) {
        throw new ReviewedStingException("containsValue() not implemented for StratificationManager");
    }

    @Override
    public V put(final List<Object> objects, final V v) {
        throw new ReviewedStingException("put() not implemented for StratificationManager");
    }

    @Override
    public V remove(final Object o) {
        throw new ReviewedStingException("remove() not implemented for StratificationManager");
    }

    @Override
    public void putAll(final Map<? extends List<Object>, ? extends V> map) {
        throw new ReviewedStingException("clear() not implemented for StratificationManager");
    }

    @Override
    public void clear() {
        throw new ReviewedStingException("clear() not implemented for StratificationManager");
    }

    @Override
    public Set<List<Object>> keySet() {
        throw new ReviewedStingException("Not yet implemented");
    }

    @Override
    public Set<Entry<List<Object>, V>> entrySet() {
        throw new ReviewedStingException("Not yet implemented");
    }

    // -------------------------------------------------------------------------------------
    //
    // utilities
    //
    // -------------------------------------------------------------------------------------

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
