/*
 * Copyright (c) 2011, The Broad Institute
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
package org.broadinstitute.sting.gatk.walkers.diffengine;


// the imports for unit testing.


import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DifferenceUnitTest extends BaseTest {
    // --------------------------------------------------------------------------------
    //
    // testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest extends TestDataProvider {
        public DiffElement tree1, tree2;
        public String difference;

        private DifferenceTest(String tree1, String tree2, String difference) {
            this(DiffNode.fromString(tree1), DiffNode.fromString(tree2), difference);
        }

        private DifferenceTest(DiffElement tree1, DiffElement tree2, String difference) {
            super(DifferenceTest.class);
            this.tree1 = tree1;
            this.tree2 = tree2;
            this.difference = difference;
        }

        public String toString() {
            return String.format("tree1=%s tree2=%s diff=%s",
                    tree1 == null ? "null" : tree1.toOneLineString(),
                    tree2 == null ? "null" : tree2.toOneLineString(),
                    difference);
        }
    }

    @DataProvider(name = "data")
    public Object[][] createTrees() {
        new DifferenceTest("A=X", "A=Y", "A:1:X!=Y");
        new DifferenceTest("A=Y", "A=X", "A:1:Y!=X");
        new DifferenceTest(DiffNode.fromString("A=X"), null, "A:1:X!=MISSING");
        new DifferenceTest(null, DiffNode.fromString("A=X"), "A:1:MISSING!=X");
        return DifferenceTest.getTests(DifferenceTest.class);
    }

    @Test(enabled = true, dataProvider = "data")
    public void testDiffToString(DifferenceTest test) {
        logger.warn("Test tree1: " + (test.tree1 == null ? "null" : test.tree1.toOneLineString()));
        logger.warn("Test tree2: " + (test.tree2 == null ? "null" : test.tree2.toOneLineString()));
        logger.warn("Test expected diff : " + test.difference);
        Difference diff = new Difference(test.tree1, test.tree2);
        logger.warn("Observed diffs     : " + diff);
        Assert.assertEquals(diff.toString(), test.difference, "Observed diff string " + diff + " not equal to expected difference string " + test.difference );

    }
}