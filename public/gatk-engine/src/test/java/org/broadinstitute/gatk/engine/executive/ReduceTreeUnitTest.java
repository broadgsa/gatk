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

package org.broadinstitute.gatk.engine.executive;


import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.Test;
import org.testng.annotations.BeforeMethod;

import org.broadinstitute.gatk.utils.BaseTest;

import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ExecutionException;
import java.util.List;
import java.util.ArrayList;
/**
 * User: hanna
 * Date: Apr 29, 2009
 * Time: 10:40:49 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Make sure the reduce tree organizes reduces in the correct way.
 */

public class ReduceTreeUnitTest extends BaseTest implements ReduceTree.TreeReduceNotifier {

    /**
     * The tree indicating reduce order.
     */
    private ReduceTree reduceTree = null;

    /**
     * 
     */
    private List<List<Integer>> reduces = new ArrayList<List<Integer>>();

    @BeforeMethod
    public void createTree() {
        reduceTree = new ReduceTree( this );
    }

    @AfterMethod
    public void destroyTree() {
        reduceTree = null;
        reduces.clear();
    }

    @Test
    public void testNoValueReduce()
        throws InterruptedException, ExecutionException {
        reduceTree.complete();
        Assert.assertEquals(reduceTree.getResult(), null, "Single-value reduce failed");
    }

    @Test
    public void testSingleValueReduce()
            throws InterruptedException, ExecutionException {
        reduceTree.addEntry( getReduceTestEntry(1) );
        reduceTree.complete();
        Assert.assertEquals(reduceTree.getResult().get(), 1, "Single-value reduce failed");
    }

    @Test(expectedExceptions=IllegalStateException.class)
    public void testIncompleteReduce()
            throws InterruptedException, ExecutionException {
        reduceTree.addEntry( getReduceTestEntry(1) );
        reduceTree.getResult().get();
    }

    @Test
    public void testDualValueReduce()
        throws InterruptedException, ExecutionException {
        reduceTree.addEntry( getReduceTestEntry(1) );
        reduceTree.addEntry( getReduceTestEntry(2) );
        reduceTree.complete();

        List<Integer> expected = new ArrayList<Integer>();
        expected.add( 1 );
        expected.add( 2 );

        // Test the result
        Assert.assertEquals(reduceTree.getResult().get(), expected, "Dual-value reduce failed");

        // Test the intermediate steps
        Assert.assertEquals(reduces.size(), 1, "Size of incoming tree reduces incorrect");
        Assert.assertEquals(reduces.get(0), expected, "Incoming tree reduce incorrect");
    }

    @Test
    public void testThreeValueReduce()
        throws InterruptedException, ExecutionException {
        List<Integer> firstExpected = new ArrayList<Integer>();
        firstExpected.add(1);
        firstExpected.add(2);

        List<Integer> finalExpected = new ArrayList<Integer>();
        finalExpected.addAll( firstExpected );
        finalExpected.add(3);

        reduceTree.addEntry( getReduceTestEntry(1) );

        Assert.assertEquals(reduces.size(), 0, "Reduce queue should be empty after entering a single element");

        reduceTree.addEntry( getReduceTestEntry(2) );

        Assert.assertEquals(reduces.size(), 1, "Reduce queue should have one element after two entries");
        Assert.assertEquals(reduces.get(0), firstExpected, "Reduce queue element is incorrect after two entries");

        reduceTree.addEntry( getReduceTestEntry(3) );

        Assert.assertEquals(reduces.size(), 1, "Reduce queue should have one element after three entries");
        Assert.assertEquals(reduces.get(0), firstExpected, "Reduce queue element is incorrect after three entries");

        reduceTree.complete();

        // Test the result
        Assert.assertEquals(reduceTree.getResult().get(), finalExpected, "Three value reduce failed");

        Assert.assertEquals(reduces.size(), 2, "Reduce queue should have two elements after three entries (complete)");
        Assert.assertEquals(reduces.get(0), firstExpected, "Reduce queue element is incorrect after three entries");
        Assert.assertEquals(reduces.get(1), finalExpected, "Reduce queue element is incorrect after three entries");
    }

    @Test
    public void testFourValueReduce()
        throws InterruptedException, ExecutionException {
        List<Integer> lhsExpected = new ArrayList<Integer>();
        lhsExpected.add(1);
        lhsExpected.add(2);

        List<Integer> rhsExpected = new ArrayList<Integer>();
        rhsExpected.add(3);
        rhsExpected.add(4);

        List<Integer> finalExpected = new ArrayList<Integer>();
        finalExpected.addAll(lhsExpected);
        finalExpected.addAll(rhsExpected);

        reduceTree.addEntry( getReduceTestEntry(1) );

        Assert.assertEquals(reduces.size(), 0, "Reduce queue should be empty after entering a single element");

        reduceTree.addEntry( getReduceTestEntry(2) );

        Assert.assertEquals(reduces.size(), 1, "Reduce queue should have one element after two entries");
        Assert.assertEquals(reduces.get(0), lhsExpected, "Reduce queue element is incorrect after two entries");

        reduceTree.addEntry( getReduceTestEntry(3) );

        Assert.assertEquals(reduces.size(), 1, "Reduce queue should have one element after three entries");
        Assert.assertEquals(reduces.get(0), lhsExpected, "Reduce queue element is incorrect after three entries");

        reduceTree.addEntry( getReduceTestEntry(4) );

        Assert.assertEquals(reduces.size(), 3, "Reduce queue should have three elements after four entries");
        Assert.assertEquals(reduces.get(0), lhsExpected, "Reduce queue element 0 is incorrect after three entries");
        Assert.assertEquals(reduces.get(1), rhsExpected, "Reduce queue element 1 is incorrect after three entries");
        Assert.assertEquals(reduces.get(2), finalExpected, "Reduce queue element 2 is incorrect after three entries");

        reduceTree.complete();

                // Test the result
        Assert.assertEquals(reduceTree.getResult().get(), finalExpected, "Four-valued reduce failed");

        // Test the working tree
        Assert.assertEquals(reduces.size(), 3, "Didn't see correct number of reduces");
        Assert.assertEquals(reduces.get(0), lhsExpected, "lhs of four value reduce failed");
        Assert.assertEquals(reduces.get(1), rhsExpected, "rhs of four value reduce failed");
        Assert.assertEquals(reduces.get(2), finalExpected, "final value four value reduce failed");
    }


    private Future getReduceTestEntry( Object value ) {
        // Create a task and run it, assuring that the tests won't block on a get.
        FutureTask task = new FutureTask( new ReduceTestEntry( value ) );
        task.run();
        return task;
    }

    public Future notifyReduce( Future lhs, Future rhs )  {
        List<Integer> reduce = new ArrayList<Integer>();

        try {
            if( lhs == null && rhs == null )
                throw new IllegalStateException("lhs and rhs are null");

            if( lhs.get() instanceof List )
                reduce.addAll((List)lhs.get());
            else
                reduce.add((Integer)lhs.get());

            if( rhs != null ) {
                if( rhs.get() instanceof List )
                    reduce.addAll((List)rhs.get());
                else
                    reduce.add((Integer)rhs.get());
            }
        }
        catch( Exception ex ) {
            // just rethrow any exceptions
            throw new RuntimeException(ex);
        }

        reduces.add( reduce );

        return getReduceTestEntry( reduce );
    }

    private class ReduceTestEntry implements Callable {
        private Object data;

        public ReduceTestEntry( Object data ) {
            this.data = data;
        }

        public Object call() {
            return data;
        }
    }
}
