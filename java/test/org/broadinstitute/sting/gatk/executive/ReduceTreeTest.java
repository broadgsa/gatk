package org.broadinstitute.sting.gatk.executive;

import org.junit.Assert;
import org.junit.Test;
import org.junit.Before;
import org.junit.After;
import org.junit.Ignore;

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

public class ReduceTreeTest implements ReduceTree.TreeReduceNotifier {

    /**
     * The tree indicating reduce order.
     */
    private ReduceTree reduceTree = null;

    /**
     * 
     */
    private List<List<Integer>> reduces = new ArrayList<List<Integer>>();

    @Before
    public void createTree() {
        reduceTree = new ReduceTree( this );
    }

    @After
    public void destroyTree() {
        reduceTree = null;
        reduces.clear();
    }

    @Test
    public void testNoValueReduce()
        throws InterruptedException, ExecutionException {
        reduceTree.complete();
        Assert.assertEquals("Single-value reduce failed", null, reduceTree.getResult());
    }

    @Test
    public void testSingleValueReduce()
            throws InterruptedException, ExecutionException {
        reduceTree.addEntry( getReduceTestEntry(1) );
        reduceTree.complete();
        Assert.assertEquals("Single-value reduce failed", 1, reduceTree.getResult().get());
    }

    @Test(expected=IllegalStateException.class)
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
        Assert.assertEquals("Dual-value reduce failed", expected, reduceTree.getResult().get());

        // Test the intermediate steps
        Assert.assertEquals("Size of incoming tree reduces incorrect", 1, reduces.size() );
        Assert.assertEquals("Incoming tree reduce incorrect", expected, reduces.get(0) );
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

        Assert.assertEquals("Reduce queue should be empty after entering a single element", 0, reduces.size());

        reduceTree.addEntry( getReduceTestEntry(2) );

        Assert.assertEquals("Reduce queue should have one element after two entries", 1, reduces.size());
        Assert.assertEquals("Reduce queue element is incorrect after two entries", firstExpected, reduces.get(0));

        reduceTree.addEntry( getReduceTestEntry(3) );

        Assert.assertEquals("Reduce queue should have one element after three entries", 1, reduces.size());
        Assert.assertEquals("Reduce queue element is incorrect after three entries", firstExpected, reduces.get(0));

        reduceTree.complete();

        // Test the result
        Assert.assertEquals("Three value reduce failed", finalExpected, reduceTree.getResult().get());

        Assert.assertEquals("Reduce queue should have two elements after three entries (complete)", 2, reduces.size());
        Assert.assertEquals("Reduce queue element is incorrect after three entries", firstExpected, reduces.get(0));
        Assert.assertEquals("Reduce queue element is incorrect after three entries", finalExpected, reduces.get(1));        
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

        Assert.assertEquals("Reduce queue should be empty after entering a single element", 0, reduces.size());

        reduceTree.addEntry( getReduceTestEntry(2) );

        Assert.assertEquals("Reduce queue should have one element after two entries", 1, reduces.size());
        Assert.assertEquals("Reduce queue element is incorrect after two entries", lhsExpected, reduces.get(0));

        reduceTree.addEntry( getReduceTestEntry(3) );

        Assert.assertEquals("Reduce queue should have one element after three entries", 1, reduces.size());
        Assert.assertEquals("Reduce queue element is incorrect after three entries", lhsExpected, reduces.get(0));

        reduceTree.addEntry( getReduceTestEntry(4) );

        Assert.assertEquals("Reduce queue should have three elements after four entries", 3, reduces.size());
        Assert.assertEquals("Reduce queue element 0 is incorrect after three entries", lhsExpected, reduces.get(0));
        Assert.assertEquals("Reduce queue element 1 is incorrect after three entries", rhsExpected, reduces.get(1));
        Assert.assertEquals("Reduce queue element 2 is incorrect after three entries", finalExpected, reduces.get(2));                

        reduceTree.complete();

                // Test the result
        Assert.assertEquals("Four-valued reduce failed",finalExpected,reduceTree.getResult().get());

        // Test the working tree
        Assert.assertEquals("Didn't see correct number of reduces", 3, reduces.size());
        Assert.assertEquals("lhs of four value reduce failed", lhsExpected, reduces.get(0));
        Assert.assertEquals("rhs of four value reduce failed", rhsExpected, reduces.get(1));
        Assert.assertEquals("final value four value reduce failed", finalExpected, reduces.get(2));
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

    @Ignore
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
