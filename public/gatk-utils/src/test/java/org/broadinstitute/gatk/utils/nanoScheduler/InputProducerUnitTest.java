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
public class InputProducerUnitTest extends BaseTest {
    @DataProvider(name = "InputProducerTest")
    public Object[][] createInputProducerTest() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( final int nElements : Arrays.asList(0, 1, 10, 100, 1000, 10000, 100000) ) {
            for ( final int queueSize : Arrays.asList(1, 10, 100) ) {
                tests.add(new Object[]{ nElements, queueSize });
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "InputProducerTest", timeOut = NanoSchedulerUnitTest.NANO_SCHEDULE_MAX_RUNTIME)
    public void testInputProducer(final int nElements, final int queueSize) throws InterruptedException {
        final List<Integer> elements = new ArrayList<Integer>(nElements);
        for ( int i = 0; i < nElements; i++ ) elements.add(i);

        final InputProducer<Integer> ip = new InputProducer<Integer>(elements.iterator());

        Assert.assertFalse(ip.allInputsHaveBeenRead(), "InputProvider said that all inputs have been read, but I haven't started reading yet");
        Assert.assertEquals(ip.getNumInputValues(), -1, "InputProvider told me that the queue was done, but I haven't started reading yet");

        int lastValue = -1;
        int nRead = 0;
        while ( ip.hasNext() ) {
            final int nTotalElements = ip.getNumInputValues();

            if ( nRead < nElements )
                Assert.assertEquals(nTotalElements, -1, "getNumInputValues should have returned -1 with not all elements read");
            // note, cannot test else case because elements input could have emptied between calls

            final InputProducer<Integer>.InputValue value = ip.next();
            if ( value.isEOFMarker() ) {
                Assert.assertEquals(nRead, nElements, "Number of input values " + nRead + " not all that are expected " + nElements);
                break;
            } else {
                Assert.assertTrue(lastValue < value.getValue(), "Read values coming out of order!");
                final int expected = lastValue + 1;
                Assert.assertEquals((int)value.getValue(), expected, "Value observed " + value.getValue() + " not equal to the expected value " + expected);
                nRead++;
                lastValue = value.getValue();
            }
        }

        Assert.assertTrue(ip.allInputsHaveBeenRead(), "InputProvider said that all inputs haven't been read, but I read them all");
        Assert.assertEquals(ip.getNumInputValues(), nElements, "Wrong number of total elements getNumInputValues");
    }
}
