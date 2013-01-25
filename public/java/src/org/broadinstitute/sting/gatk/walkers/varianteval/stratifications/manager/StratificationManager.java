/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Represents the full state space of all stratification combinations
 *
 * @author Mark DePristo
 * @since 3/27/12
 */
public class StratificationManager<K extends Stratifier, V> implements Map<List<Object>, V> {
    private final StratNode<K> root;
    private final int size;

    private final ArrayList<K> stratifiers;

    // values associated with each key
    private final ArrayList<V> valuesByKey;
    private final ArrayList<List<Object>> stratifierValuesByKey;
    private final ArrayList<String> keyStrings;

    // -------------------------------------------------------------------------------------
    //
    // creating the manager
    //
    // -------------------------------------------------------------------------------------

    /**
     * Create a new StratificationManager with nodes to store data for all combinations
     * of the ordered list of strats
     *
     * @param strats ordered list of stratifications to representation
     */
    @Requires("!strats.isEmpty()")
    public StratificationManager(final List<K> strats) {
        this.stratifiers = new ArrayList<K>(strats);

        // construct and store the full tree of strats
        this.root = buildStratificationTree(new LinkedList<K>(strats));
        // assign the linear key ordering to the leafs
        assignKeys(root);

        // cache the size, and check for a bad state
        this.size = root.size();
        if ( this.size == 0 )
            throw new ReviewedStingException("Size == 0 in StratificationManager");

        // prepare the assocated data vectors mapping from key -> data
        this.valuesByKey = new ArrayList<V>(size());
        this.stratifierValuesByKey = new ArrayList<List<Object>>(size());
        this.keyStrings = new ArrayList<String>(size());
        for ( int i = 0; i < size(); i++ ) {
            this.valuesByKey.add(null);
            this.stratifierValuesByKey.add(null);
            this.keyStrings.add(null);
        }

        assignStratifierValuesByKey(root);
    }

    /**
     * Recursive construction helper for main constructor that fills into the
     * complete tree of StratNodes.  This function returns the complete tree
     * suitable for associating data with each combinatino of keys.  Note
     * that the tree is not fully complete as the keys are not yet set for
     * each note (see assignStratifierValuesByKey)
     *
     * @param strats
     * @return
     */
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

    /**
     * Set the key for each leaf from root, in order from 0 to N - 1 for N leaves in the tree
     * @param root
     */
    @Requires("root == this.root")
    private void assignKeys(final StratNode<K> root) {
        int key = 0;
        for ( final StratNode<K> node : root ) {
            if ( node.isLeaf() )
                node.setKey(key++);
        }
    }

    /**
     * Entry point to recursive tool that fills in the list of state values corresponding
     * to each key.  After this function is called you can map from key -> List of StateValues
     * instead of walking the tree to find the key and reading the list of state values
     *
     * @param root
     */
    private void assignStratifierValuesByKey(final StratNode<K> root) {
        assignStratifierValuesByKey(root, new LinkedList<Object>());

        // do a last sanity check that no key has null value after assigning
        for ( List<Object> stateValues : stratifierValuesByKey )
            if ( stateValues == null )
                throw new ReviewedStingException("Found a null state value set that's null");
    }

    private void assignStratifierValuesByKey(final StratNode<K> node, final LinkedList<Object> states) {
        if ( node.isLeaf() ) { // we're here!
            if ( states.isEmpty() )
                throw new ReviewedStingException("Found a leaf node with an empty state values vector");
            stratifierValuesByKey.set(node.getKey(), Collections.unmodifiableList(new ArrayList<Object>(states)));
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

    /**
     * How many states are held in this stratification manager?
     * @return
     */
    @Ensures("result >= 0")
    public int size() {
        return size;
    }

    @Ensures("result != null")
    protected StratNode<K> getRoot() {
        return root;
    }

    @Ensures("result != null")
    public List<K> getStratifiers() {
        return stratifiers;
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

    public List<Object> getStatesForKey(final int key) {
        final List<Object> states = new ArrayList<Object>(stratifiers.size());
        for ( int i = 0; i < stratifiers.size(); i++ ) {
            final Object stratValue = stratifierValuesByKey.get(key).get(i);
            states.add(stratValue);
        }
        return states;
    }

    public List<Pair<K, Object>> getStratsAndStatesForKey(final int key) {
        final List<Pair<K, Object>> states = new ArrayList<Pair<K, Object>>(stratifiers.size());
        for ( int i = 0; i < stratifiers.size(); i++ ) {
            final K strat = stratifiers.get(i);
            final Object stratValue = stratifierValuesByKey.get(key).get(i);
            states.add(new Pair<K, Object>(strat, stratValue));
        }
        return states;
    }

    public String getStratsAndStatesStringForKey(final int key) {
        if ( keyStrings.get(key) == null ) {
            StringBuilder b = new StringBuilder();
            for ( int i = 0; i < stratifiers.size(); i++ ) {
                final K strat = stratifiers.get(i);
                final Object stratValue = stratifierValuesByKey.get(key).get(i);
                b.append(strat.toString()).append(":").append(stratValue.toString());
            }
            keyStrings.set(key, b.toString());
        }
        
        return keyStrings.get(key);
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
        final List<List<Object>> combined = new ArrayList<List<Object>>(first.size());
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

    public interface Combiner<V> {
        /** take two values of type V and return a combined value of type V */
        public V combine(final V lhs, final V rhs);
    }

    /**
     * Remaps the stratifications from one stratification set to another, combining
     * the values in V according to the combiner function.
     *
     * stratifierToReplace defines a set of states S1, while newStratifier defines
     * a new set S2.  remappedStates is a map from all of S1 into at least some of
     * S2.  This function creates a new, fully initialized manager where all of the
     * data in this new manager is derived from the original data in this object
     * combined according to the mapping remappedStates.  When multiple
     * elements of S1 can map to the same value in S2, these are sequentially
     * combined by the function combiner.  Suppose for example at states s1, s2, and
     * s3 all map to N1.  Eventually the value associated with state N1 would be
     *
     *   value(N1) = combine(value(s1), combine(value(s2), value(s3))
     *
     * in some order for s1, s2, and s3, which is not defined.  Note that this function
     * only supports combining one stratification at a time, but in principle a loop over
     * stratifications and this function could do the multi-dimensional collapse.
     *
     * @param stratifierToReplace
     * @param newStratifier
     * @param combiner
     * @param remappedStates
     * @return
     */
    public StratificationManager<K, V> combineStrats(final K stratifierToReplace,
                                                     final K newStratifier,
                                                     final Combiner<V> combiner,
                                                     final Map<Object, Object> remappedStates) {
        // make sure the mapping is reasonable
        if ( ! newStratifier.getAllStates().containsAll(remappedStates.values()) )
            throw new ReviewedStingException("combineStrats: remapped states contains states not found in newStratifer state set");

        if ( ! remappedStates.keySet().containsAll(stratifierToReplace.getAllStates()) )
            throw new ReviewedStingException("combineStrats: remapped states missing mapping for some states");

        // the new strats are the old ones with the single replacement
        final List<K> newStrats = new ArrayList<K>(getStratifiers());
        final int stratOffset = newStrats.indexOf(stratifierToReplace);
        if ( stratOffset == -1 )
            throw new ReviewedStingException("Could not find strat to replace " + stratifierToReplace + " in existing strats " + newStrats);
        newStrats.set(stratOffset, newStratifier);

        // create an empty but fully initialized new manager
        final StratificationManager<K, V> combined = new StratificationManager<K, V>(newStrats);

        // for each key, get its state, update it according to the map, and update the combined manager
        for ( int key = 0; key < size(); key++ ) {
            // the new state is just the old one with the replacement
            final List<Object> newStates = new ArrayList<Object>(getStatesForKey(key));
            final Object oldState = newStates.get(stratOffset);
            final Object newState = remappedStates.get(oldState);
            newStates.set(stratOffset, newState);

            // look up the new key given the new state
            final int combinedKey = combined.getKey(newStates);
            if ( combinedKey == -1 ) throw new ReviewedStingException("Couldn't find key for states: " + Utils.join(",", newStates));

            // combine the old value with whatever new value is in combined already
            final V combinedValue = combiner.combine(combined.get(combinedKey), get(key));

            // update the value associated with combined key
            combined.set(combinedKey, combinedValue);
        }

        return combined;
    }
}