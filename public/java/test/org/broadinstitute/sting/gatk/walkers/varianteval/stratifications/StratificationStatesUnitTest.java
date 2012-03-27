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

// our package
package org.broadinstitute.sting.gatk.walkers.varianteval.stratifications;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.FileNotFoundException;
import java.util.*;


public class StratificationStatesUnitTest extends BaseTest {
    @BeforeClass
    public void init() throws FileNotFoundException {
    }

    // --------------------------------------------------------------------------------
    //
    // Basic tests Provider
    //
    // --------------------------------------------------------------------------------

    private class StratificationStatesTestProvider extends TestDataProvider {
        final List<List<Integer>> allStates;
        final List<ListAsSetOfStates> asSetOfStates = new ArrayList<ListAsSetOfStates>();
        final int nStates;
        
        public StratificationStatesTestProvider(final List<Integer> ... allStates) {
            super(StratificationStatesTestProvider.class);
            this.allStates = Arrays.asList(allStates);

            int nStates = 1;
            for ( List<Integer> states : this.allStates ) { 
                nStates *= states.size();
                asSetOfStates.add(new ListAsSetOfStates(states));
            }
            this.nStates = nStates;
        }
//        private String getName() {
//            return String.format("probs=%s expectedRegions=%s", Utils.join(",", probs), Utils.join(",", expectedRegions));
//        }

        public List<ListAsSetOfStates> getStateSpaceList() {
            return asSetOfStates;
        }
        
        public Queue<List<String>> getAllCombinations() {
            return getAllCombinations(new LinkedList<List<Integer>>(allStates));
        }

        private Queue<List<String>> getAllCombinations(Queue<List<Integer>> states) {
            if ( states.isEmpty() ) 
                return new LinkedList<List<String>>();
            else {
                List<Integer> head = states.poll();
                Queue<List<String>> substates = getAllCombinations(states);
                Queue<List<String>> newStates = new LinkedList<List<String>>();
                for ( int e : head) {
                    for ( List<String> state : substates ) {
                        List<String> newState = new LinkedList<String>();
                        newState.add(Integer.toString(e));
                        newState.addAll(state);
                        newStates.add(newState);
                    }
                }
                return newStates;
            }
        }
    }

    private class ListAsSetOfStates implements SetOfStates {
        final List<String> integers;

        private ListAsSetOfStates(final List<Integer> integers) {
            this.integers = new ArrayList<String>(integers.size());
            for ( int i : integers )
                this.integers.add(Integer.toString(i));
        }

        @Override
        public List<String> getAllStates() {
            return integers;
        }
    }

    @DataProvider(name = "StratificationStatesTestProvider")
    public Object[][] makeStratificationStatesTestProvider() {
        new StratificationStatesTestProvider(Arrays.asList(0));
        new StratificationStatesTestProvider(Arrays.asList(0, 1));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4), Arrays.asList(5, 6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3, 4, 5), Arrays.asList(6, 7));
        new StratificationStatesTestProvider(Arrays.asList(0, 1), Arrays.asList(2, 3), Arrays.asList(4, 5), Arrays.asList(6, 7));
        return StratificationStatesTestProvider.getTests(StratificationStatesTestProvider.class);
    }

    @Test(dataProvider = "StratificationStatesTestProvider")
    public void testStratificationStatesTestProvider(StratificationStatesTestProvider cfg) {
        StratificationStates<ListAsSetOfStates> stratificationStates = new StratificationStates<ListAsSetOfStates>(cfg.getStateSpaceList());

        Assert.assertEquals(stratificationStates.getNStates(), cfg.nStates);
        
        int nLeafs = 0;
        for ( final StratNode node : stratificationStates.getRoot() ) {
            if ( node.isLeaf() )
                nLeafs++;
        }
        Assert.assertEquals(nLeafs, cfg.nStates, "Unexpected number of leaves");
        
        Set<Integer> seenKeys = new HashSet<Integer>(cfg.nStates);
        for ( final StratNode node : stratificationStates.getRoot() ) {
            if ( node.isLeaf() ) {
                Assert.assertFalse(seenKeys.contains(node.getKey()), "Already seen the key");
                seenKeys.add(node.getKey());
            }
        }

        seenKeys.clear();
        for ( List<String> state : cfg.getAllCombinations() ) {
            final int key = stratificationStates.getKey(state);
            Assert.assertFalse(seenKeys.contains(key), "Already saw state mapping to this key");
            seenKeys.add(key);
        }
    }
}