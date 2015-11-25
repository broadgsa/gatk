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

package org.broadinstitute.gatk.engine.datasources.providers;

import org.testng.Assert;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.annotations.BeforeMethod;


import org.testng.annotations.Test;
import org.broadinstitute.gatk.utils.BaseTest;

import java.util.Collection;
import java.util.Collections;
import java.util.Arrays;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 1:56:02 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test basic functionality of the shard data provider.
 */

public class ShardDataProviderUnitTest extends BaseTest {
    /**
     * Provider to test.  Should be recreated for every test.
     */
    private ShardDataProvider provider = null;

    @BeforeMethod
    public void createProvider() {
        provider = new LocusShardDataProvider( null,null,null,null,null,null,null );
    }

    /**
     * Test whether views are closed when the provider closes.
     */
    @Test
    public void testClose() {
        TestView testView = new TestView( provider );
        Assert.assertFalse(testView.closed,"View is currently closed but should be open");

        provider.close();
        Assert.assertTrue(testView.closed,"View is currently open but should be closed");
    }

    /**
     * Test whether multiple of the same view can be registered and all get a close method.
     */
    @Test
    public void testMultipleClose() {
        Collection<TestView> testViews = Arrays.asList(new TestView(provider),new TestView(provider));
        for( TestView testView: testViews )
            Assert.assertFalse(testView.closed,"View is currently closed but should be open");

        provider.close();
        for( TestView testView: testViews )
            Assert.assertTrue(testView.closed,"View is currently open but should be closed");
    }

    /**
     * Try adding a view which conflicts with some other view that's already been registered.
     */
    @Test(expectedExceptions= ReviewedGATKException.class)
    public void testAddViewWithExistingConflict() {
        View initial = new ConflictingTestView( provider );
        View conflictsWithInitial = new TestView( provider );
    }

    /**
     * Try adding a view which has a conflict with a previously registered view.
     */
    @Test(expectedExceptions= ReviewedGATKException.class)
    public void testAddViewWithNewConflict() {
        View conflictsWithInitial = new TestView( provider );
        View initial = new ConflictingTestView( provider );
    }

    /**
     * A simple view for testing interactions between views attached to the ShardDataProvider.
     */
    private class TestView implements View {
        /**
         * Is the test view currently closed.
         */
        private boolean closed = false;

        /**
         * Create a new test view wrapping the given provider.
         * @param provider
         */
        public TestView( ShardDataProvider provider ) {
            provider.register(this);            
        }

        /**
         * Gets conflicting views.  In this case, none conflict.
         * @return
         */
        public Collection<Class<? extends View>> getConflictingViews() { return Collections.emptyList(); }

        /**
         * Close this view.
         */
        public void close() { this.closed = true; }
    }

    /**
     * Another view that conflicts with the one above.
     */
    private class ConflictingTestView implements View {
        public ConflictingTestView( ShardDataProvider provider ) { provider.register(this); }

        public Collection<Class<? extends View>> getConflictingViews() {
            return Collections.<Class<? extends View>>singleton(TestView.class);
        }

        public void close() {}
    }
}
