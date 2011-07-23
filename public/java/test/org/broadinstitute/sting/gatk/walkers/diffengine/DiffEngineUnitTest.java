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

import java.util.*;

/**
 * Basic unit test for DifferableReaders in reduced reads
 */
public class DiffEngineUnitTest extends BaseTest {
    DiffEngine engine;

    @BeforeClass(enabled = true)
    public void createDiffEngine() {
        engine = new DiffEngine();
    }

    // --------------------------------------------------------------------------------
    //
    // Difference testing routines
    //
    // --------------------------------------------------------------------------------

    private class DifferenceTest extends TestDataProvider {
        public DiffElement tree1, tree2;
        public List<String> differences;

        private DifferenceTest(String tree1, String tree2) {
            this(tree1, tree2, Collections.<String>emptyList());
        }

        private DifferenceTest(String tree1, String tree2, String difference) {
            this(tree1, tree2, Arrays.asList(difference));
        }

        private DifferenceTest(String tree1, String tree2, List<String> differences) {
            super(DifferenceTest.class);
            this.tree1 = DiffNode.fromString(tree1);
            this.tree2 = DiffNode.fromString(tree2);
            this.differences = differences;
        }

        public String toString() {
            return String.format("tree1=%s tree2=%s diff=%s",
                    tree1.toOneLineString(), tree2.toOneLineString(), differences);
        }
    }

    @DataProvider(name = "trees")
    public Object[][] createTrees() {
        new DifferenceTest("A=X", "A=X");
        new DifferenceTest("A=X", "A=Y", "A:X!=Y");
        new DifferenceTest("A=X", "B=X", Arrays.asList("A:X!=MISSING", "B:MISSING!=X"));
        new DifferenceTest("A=(X=1)", "B=(X=1)", Arrays.asList("A:(X=1)!=MISSING", "B:MISSING!=(X=1)"));
        new DifferenceTest("A=(X=1)", "A=(X=1)");
        new DifferenceTest("A=(X=1 Y=2)", "A=(X=1 Y=2)");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=3))");
        new DifferenceTest("A=(X=1)", "A=(X=2)", "A.X:1!=2");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2 B=(Z=4))", "A.B.Z:3!=4");
        new DifferenceTest("A=(X=1)", "A=(X=1 Y=2)", "A.Y:MISSING!=2");
        new DifferenceTest("A=(X=1 Y=2 B=(Z=3))", "A=(X=1 Y=2)", "A.B:(Z=3)!=MISSING");
        return DifferenceTest.getTests(DifferenceTest.class);
    }

    @Test(enabled = true, dataProvider = "trees")
    public void testDiffs(DifferenceTest test) {
        logger.warn("Test tree1: " + test.tree1.toOneLineString());
        logger.warn("Test tree2: " + test.tree2.toOneLineString());

        List<Difference> diffs = engine.diff(test.tree1, test.tree2);
        logger.warn("Test expected diff : " + test.differences);
        logger.warn("Observed diffs     : " + diffs);
    }

    // --------------------------------------------------------------------------------
    //
    // Low-level routines for summarizing differences
    //
    // --------------------------------------------------------------------------------

    @Test(enabled = true)
    public void testLongestCommonPostfix() {
        testLongestCommonPostfixHelper("A", "A", 1);
        testLongestCommonPostfixHelper("A", "B", 0);
        testLongestCommonPostfixHelper("A.B", "A.B", 2);
        testLongestCommonPostfixHelper("A.B.C", "A.B.C", 3);
        testLongestCommonPostfixHelper("A.B.C", "X.B.C", 2);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "X.Y.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.X.C", 1);
        testLongestCommonPostfixHelper("A.B.C", "A.X.Z", 0);
        testLongestCommonPostfixHelper("A.B.C", "A.B.Z", 0);
    }

    public void testLongestCommonPostfixHelper(String p1, String p2, int expected) {
        String[] parts1 = p1.split("\\.");
        String[] parts2 = p2.split("\\.");
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        Assert.assertEquals(obs, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }

    @Test(enabled = true, dependsOnMethods = "testLongestCommonPostfix")
    public void testSummarizePath() {
        testSummarizePathHelper("A", "A", "A");
        testSummarizePathHelper("A", "B", "*");
        testSummarizePathHelper("A.B", "A.B", "A.B");
        testSummarizePathHelper("A.B", "X.B", "*.B");
        testSummarizePathHelper("A.B", "X.Y", "*.*");
        testSummarizePathHelper("A.B.C", "A.B.C", "A.B.C");
        testSummarizePathHelper("A.B.C", "X.B.C", "*.B.C");
        testSummarizePathHelper("A.B.C", "X.Y.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "X.Y.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.X.C", "*.*.C");
        testSummarizePathHelper("A.B.C", "A.X.Z", "*.*.*");
        testSummarizePathHelper("A.B.C", "A.B.Z", "*.*.*");
    }

    public void testSummarizePathHelper(String p1, String p2, String expected) {
        String[] parts1 = DiffEngine.diffNameToPath(p1);
        String[] parts2 = DiffEngine.diffNameToPath(p2);
        int obs = DiffEngine.longestCommonPostfix(parts1, parts2);
        String path = DiffEngine.summarizedPath(parts2, obs);
        Assert.assertEquals(path, expected, "p1=" + p1 + " p2=" + p2 + " failed");
    }

    // --------------------------------------------------------------------------------
    //
    // High-level difference summary
    //
    // --------------------------------------------------------------------------------

    private class SummarizeDifferenceTest extends TestDataProvider {
        List<String> diffs = new ArrayList<String>();
        List<String> expecteds = new ArrayList<String>();

        public SummarizeDifferenceTest() { super(SummarizeDifferenceTest.class); }

        public SummarizeDifferenceTest addDiff(String... diffsToAdd) {
            diffs.addAll(Arrays.asList(diffsToAdd));
            return this;
        }

        public SummarizeDifferenceTest addSummary(String... expectedSummary) {
            expecteds.addAll(Arrays.asList(expectedSummary));
            return this;
        }

        public String toString() {
            return String.format("diffs=%s => expected=%s", diffs, expecteds);
        }

        public void test() {
            List<String[]> diffPaths = new ArrayList<String[]>(diffs.size());
            for ( String diff : diffs ) { diffPaths.add(DiffEngine.diffNameToPath(diff)); }

            List<Difference> sumDiffs = engine.summarizedDifferencesOfPathsFromString(diffs);

            Assert.assertEquals(sumDiffs.size(), expecteds.size(), "Unexpected number of summarized differences: " + sumDiffs);

            for ( int i = 0; i < sumDiffs.size(); i++ ) {
                Difference sumDiff = sumDiffs.get(i);
                String expected = expecteds.get(i);
                String[] pathCount = expected.split(":");
                String path = pathCount[0];
                int count = Integer.valueOf(pathCount[1]);
                Assert.assertEquals(sumDiff.getPath(), path, "Unexpected path at: " + expected + " obs=" + sumDiff + " all=" + sumDiffs);
                Assert.assertEquals(sumDiff.getCount(), count, "Unexpected counts at: " + expected + " obs=" + sumDiff + " all=" + sumDiffs);
            }
        }
    }

    @DataProvider(name = "summaries")
    public Object[][] createSummaries() {
        new SummarizeDifferenceTest().addDiff("A", "A").addSummary("A:2");
        new SummarizeDifferenceTest().addDiff("A", "B").addSummary("A:1", "B:1");
        new SummarizeDifferenceTest().addDiff("A", "A", "A").addSummary("A:3");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B").addSummary("A:3", "B:1");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B", "B").addSummary("A:3", "B:2");
        new SummarizeDifferenceTest().addDiff("A", "A", "A", "B", "B", "C").addSummary("A:3", "B:2", "C:1");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X").addSummary("A.X:2");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X", "B.X").addSummary("*.X:3", "A.X:2", "B.X:1");
        new SummarizeDifferenceTest().addDiff("A.X", "A.X", "B.X", "B.X").addSummary("*.X:4", "A.X:2", "B.X:2");
        new SummarizeDifferenceTest().addDiff("A.B.C", "X.B.C").addSummary("*.B.C:2", "A.B.C:1", "X.B.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "X.Y.C", "X.Y.C").addSummary("*.*.C:3", "X.Y.C:2", "A.B.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "X.Y.C").addSummary("*.*.C:3", "A.B.C:1", "A.X.C:1", "X.Y.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "B.X.C").addSummary("*.*.C:3", "*.X.C:2", "A.B.C:1", "A.X.C:1", "B.X.C:1");
        new SummarizeDifferenceTest().addDiff("A.B.C", "A.X.C", "B.X.C", "B.X.C").addSummary("*.*.C:4", "*.X.C:3", "B.X.C:2", "A.B.C:1", "A.X.C:1");

        return SummarizeDifferenceTest.getTests(SummarizeDifferenceTest.class);
    }


    @Test(enabled = true, dependsOnMethods = "testSummarizePath", dataProvider = "summaries")
    public void testSummarizeDifferences(SummarizeDifferenceTest test) {
        test.test();
    }
}