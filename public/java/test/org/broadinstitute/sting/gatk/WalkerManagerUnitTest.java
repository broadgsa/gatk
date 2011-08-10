/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk;

import org.testng.Assert;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.qc.CountLociWalker;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * Tests basic functionality of the walker manager.
 */
public class WalkerManagerUnitTest {
    private static WalkerManager walkerManager;

    @BeforeClass
    public void setUp() {
        walkerManager = new WalkerManager();
    }

    @Test
    public void testPresentWalker() {
        Walker countLociWalker = walkerManager.createByName("CountLoci");
        Assert.assertEquals(CountLociWalker.class,countLociWalker.getClass());
    }

    @Test(expectedExceptions=UserException.class)
    public void testAbsentWalker() {
        walkerManager.createByName("Missing");
    }

    @Test(expectedExceptions=DynamicClassResolutionException.class)
    public void testUninstantiableWalker() {
        walkerManager.createByName("Uninstantiable");
    }
}

@Hidden
class UninstantiableWalker extends Walker<Integer,Long> {
    // Private constructor will generate uninstantiable message
    private UninstantiableWalker() {}
    public Long reduceInit() { return 0L; }
    public Long reduce(Integer value, Long accum) { return 0L; }
}
