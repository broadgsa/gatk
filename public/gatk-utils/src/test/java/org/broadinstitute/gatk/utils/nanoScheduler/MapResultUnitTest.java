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

package org.broadinstitute.gatk.utils.nanoScheduler;

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
* UnitTests for the InputProducer
*
* User: depristo
* Date: 8/24/12
* Time: 11:25 AM
* To change this template use File | Settings | File Templates.
*/
public class MapResultUnitTest {
    @DataProvider(name = "CompareTester")
    public Object[][] createCompareTester() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( int id1 = 0; id1 < 10; id1++ ) {
            for ( int id2 = 0; id2 < 10; id2++ ) {
                tests.add(new Object[]{ id1, id2, Integer.valueOf(id1).compareTo(id2)});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "CompareTester")
    public void testInputProducer(final int id1, final int id2, final int comp ) throws InterruptedException {
        final MapResult<Integer> mr1 = new MapResult<Integer>(id1, id1);
        final MapResult<Integer> mr2 = new MapResult<Integer>(id2, id2);
        Assert.assertEquals(mr1.compareTo(mr2), comp, "Compare MapResultsUnitTest failed");
    }
}
