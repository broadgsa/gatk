/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.annotator;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by gauthier on 11/1/16.
 */
public class HistogramUnitTest {
    private final Double EPSILON = 0.001;

    @Test
    public void testAdd() throws Exception {
        Histogram bimodalHist = new Histogram();
        for(int i=0; i<=100; i++) {
            bimodalHist.add(1+i/1000.0);
        }
        Assert.assertEquals(bimodalHist.get(1.0), new Integer(100), "");
        Assert.assertEquals(bimodalHist.get(1.1), new Integer(1), "");
    }

    @Test
    public void testAddingQuantizedValues() throws Exception {
        Histogram hist = new Histogram();
        for(int i=0; i<100; i++) {
            hist.add(1.2);
        }
        Assert.assertEquals(hist.get(1.2), new Integer(100));
        Assert.assertEquals(hist.median(), 1.2, EPSILON);
    }

    @Test
    public void testBulkAdd() throws Exception {
        Histogram bimodalHist = new Histogram();
        for(int i=0; i<=100; i++) {
            bimodalHist.add(1+i/1000.0, 2);
        }
        Assert.assertEquals(bimodalHist.get(1.0), new Integer(200), "");
        Assert.assertEquals(bimodalHist.get(1.1), new Integer(2), "");
    }

    @Test
    public void testMedianOfEvens() throws Exception {
        Histogram bimodalHist = new Histogram();
        for(int i = 0; i<10; i++) {
            bimodalHist.add(10.0);
            bimodalHist.add(20.0);
        }

        Assert.assertEquals(bimodalHist.median(), 15.0, EPSILON, "");
    }

    @Test
    public void testMedianOfOdds() throws Exception {
        Histogram bimodalHist = new Histogram();
        for(int i = 0; i<10; i++) {
            bimodalHist.add(10.0);
            bimodalHist.add(20.0);
        }
        bimodalHist.add(20.0);

        Assert.assertEquals(bimodalHist.median(), 20.0, EPSILON, "");
    }

    @Test
    public void testMedianOfEmptyHist() throws Exception {
        Histogram empty = new Histogram();
        Assert.assertNull(empty.median());
    }

    @Test
    public void testMedianOfSingleItem() throws Exception {
        Histogram singleItem = new Histogram();
        singleItem.add(20.0);
        Assert.assertEquals(singleItem.median(), 20.0, EPSILON, "");
    }

}