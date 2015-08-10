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

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CountReadsUnitTest extends BaseTest {

    @Test
    public void testReadsDoNotOverflowInt() {

        final CountReads walker = new CountReads();

        final long moreThanMaxInt = ((long)Integer.MAX_VALUE) + 1L;

        Long sum = walker.reduceInit();

        for ( long i = 0L; i < moreThanMaxInt; i++ ) {
            final Integer x = walker.map(null, null, null);
            sum = walker.reduce(x, sum);
        }

        Assert.assertEquals(sum.longValue(), moreThanMaxInt);
        Assert.assertTrue(sum.longValue() > (long) Integer.MAX_VALUE);
    }
}
