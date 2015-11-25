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

package org.broadinstitute.gatk.utils.jna.clibrary;

import com.sun.jna.NativeLong;
import com.sun.jna.ptr.NativeLongByReference;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class LibCUnitTest extends BaseTest {

    @Test
    public void testEnvironment() {
        String testProperty = "test_property";
        String testValue = "value";
        Assert.assertEquals(LibC.getenv(testProperty), null);
        Assert.assertEquals(LibC.setenv(testProperty, testValue, 1), 0);
        Assert.assertEquals(LibC.getenv(testProperty), testValue);
        Assert.assertEquals(LibC.unsetenv(testProperty), 0);
        Assert.assertEquals(LibC.getenv(testProperty), null);
    }

    @Test
    public void testDifftime() throws Exception {
        // Pointer to hold the times
        NativeLongByReference ref = new NativeLongByReference();

        // time() returns -1 on error.
        NativeLong err = new NativeLong(-1L);

        LibC.time(ref);
        NativeLong time0 = ref.getValue();
        Assert.assertNotSame(time0, err, "Time 0 returned an error (-1).");

        Thread.sleep(5000L);

        LibC.time(ref);
        NativeLong time1 = ref.getValue();
        Assert.assertNotSame(time1, err, "Time 1 returned an error (-1).");

        Assert.assertNotSame(time1, time0, "Time 1 returned same time as Time 0.");

        double diff = LibC.difftime(time1, time0);
        Assert.assertTrue(diff >= 5, "Time difference was not greater than 5 seconds: " + diff);
    }
}
