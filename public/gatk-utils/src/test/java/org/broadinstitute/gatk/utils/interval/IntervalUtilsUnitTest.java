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

package org.broadinstitute.gatk.utils.interval;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.io.FileUtils;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.IntervalBinding;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * test out the interval utility methods
 */
public class IntervalUtilsUnitTest extends BaseTest {
    // used to seed the genome loc parser with a sequence dictionary
    private SAMFileHeader hg18Header;
    private GenomeLocParser hg18GenomeLocParser;
    private List<GenomeLoc> hg18ReferenceLocs;
    private SAMFileHeader hg19Header;
    private GenomeLocParser hg19GenomeLocParser;
    private List<GenomeLoc> hg19ReferenceLocs;
    private List<GenomeLoc> hg19exomeIntervals;

    private List<GenomeLoc> getLocs(String... intervals) {
        return getLocs(Arrays.asList(intervals));
    }

    private List<GenomeLoc> getLocs(List<String> intervals) {
        if (intervals.size() == 0)
            return hg18ReferenceLocs;
        List<GenomeLoc> locs = new ArrayList<GenomeLoc>();
        for (String interval: intervals)
            locs.add(hg18GenomeLocParser.parseGenomeLoc(interval));
        return Collections.unmodifiableList(locs);
    }

    @BeforeClass
    public void init() {
        File hg18Ref = new File(BaseTest.hg18Reference);
        try {
            final ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(hg18Ref);
            hg18Header = new SAMFileHeader();
            hg18Header.setSequenceDictionary(seq.getSequenceDictionary());
            hg18GenomeLocParser = new GenomeLocParser(seq);
            hg18ReferenceLocs = Collections.unmodifiableList(GenomeLocSortedSet.createSetFromSequenceDictionary(seq.getSequenceDictionary()).toList()) ;
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(hg18Ref,ex);
        }

        File hg19Ref = new File(BaseTest.hg19Reference);
        try {
            final ReferenceSequenceFile seq = new CachingIndexedFastaSequenceFile(hg19Ref);
            hg19Header = new SAMFileHeader();
            hg19Header.setSequenceDictionary(seq.getSequenceDictionary());
            hg19GenomeLocParser = new GenomeLocParser(seq);
            hg19ReferenceLocs = Collections.unmodifiableList(GenomeLocSortedSet.createSetFromSequenceDictionary(seq.getSequenceDictionary()).toList()) ;

            hg19exomeIntervals = Collections.unmodifiableList(IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(hg19Intervals)));
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(hg19Ref,ex);
        }
    }

    // -------------------------------------------------------------------------------------
    //
    // tests to ensure the quality of the interval cuts of the interval cutting functions
    //
    // -------------------------------------------------------------------------------------

    private class IntervalSlicingTest extends TestDataProvider {
        public int parts;
        public double maxAllowableVariance;

        private IntervalSlicingTest(final int parts, final double maxAllowableVariance) {
            super(IntervalSlicingTest.class);
            this.parts = parts;
            this.maxAllowableVariance = maxAllowableVariance;
        }

        public String toString() {
            return String.format("IntervalSlicingTest parts=%d maxVar=%.2f", parts, maxAllowableVariance);
        }
    }

    @DataProvider(name = "intervalslicingdata")
    public Object[][] createTrees() {
        new IntervalSlicingTest(1, 0);
        new IntervalSlicingTest(2, 1);
        new IntervalSlicingTest(5, 1);
        new IntervalSlicingTest(10, 1);
        new IntervalSlicingTest(67, 1);
        new IntervalSlicingTest(100, 1);
        new IntervalSlicingTest(500, 1);
        new IntervalSlicingTest(1000, 1);
        return IntervalSlicingTest.getTests(IntervalSlicingTest.class);
    }

    @Test(enabled = true, dataProvider = "intervalslicingdata")
    public void testFixedScatterIntervalsAlgorithm(IntervalSlicingTest test) {
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19exomeIntervals, test.parts);

        long totalSize = IntervalUtils.intervalSize(hg19exomeIntervals);
        long idealSplitSize = totalSize / test.parts;

        long sumOfSplitSizes = 0;
        int counter = 0;
        for ( final List<GenomeLoc> split : splits ) {
            long splitSize = IntervalUtils.intervalSize(split);
            double sigma = (splitSize - idealSplitSize) / (1.0 * idealSplitSize);
            //logger.warn(String.format("Split %d size %d ideal %d sigma %.2f", counter, splitSize, idealSplitSize, sigma));
            counter++;
            sumOfSplitSizes += splitSize;
            Assert.assertTrue(Math.abs(sigma) <= test.maxAllowableVariance, String.format("Interval %d (size %d ideal %d) has a variance %.2f outside of the tolerated range %.2f", counter, splitSize, idealSplitSize, sigma, test.maxAllowableVariance));
        }

        Assert.assertEquals(totalSize, sumOfSplitSizes, "Split intervals don't contain the exact number of bases in the origianl intervals");
    }

    // -------------------------------------------------------------------------------------
    //
    // splitLocusIntervals tests
    //
    // -------------------------------------------------------------------------------------

    /** large scale tests for many intervals */
    private class SplitLocusIntervalsTest extends TestDataProvider {
        final List<GenomeLoc> originalIntervals;
        final public int parts;

        private SplitLocusIntervalsTest(final String name, List<GenomeLoc> originalIntervals, final int parts) {
            super(SplitLocusIntervalsTest.class, name);
            this.parts = parts;
            this.originalIntervals = originalIntervals;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "IntervalRepartitionTest")
    public Object[][] createIntervalRepartitionTest() {
        for ( int parts : Arrays.asList(1, 2, 3, 10, 13, 100, 151, 1000, 10000) ) {
        //for ( int parts : Arrays.asList(10) ) {
            new SplitLocusIntervalsTest("hg19RefLocs", hg19ReferenceLocs, parts);
            new SplitLocusIntervalsTest("hg19ExomeLocs", hg19exomeIntervals, parts);
        }

        return SplitLocusIntervalsTest.getTests(SplitLocusIntervalsTest.class);
    }

    @Test(enabled = true, dataProvider = "IntervalRepartitionTest")
    public void testIntervalRepartition(SplitLocusIntervalsTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.originalIntervals, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.parts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test overall size
        final long originalSize = IntervalUtils.intervalSize(test.originalIntervals);
        final long flatSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(flatSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        // test size of each split
        final long ideal = (long)Math.floor(originalSize / (1.0 * test.parts));
        final long maxSize = ideal + (originalSize % test.parts) * test.parts; // no more than N * rounding error in size
        for ( final List<GenomeLoc> split : splitByLocus ) {
            final long splitSize = IntervalUtils.intervalSize(split);
            Assert.assertTrue(splitSize >= ideal && splitSize <= maxSize,
                    String.format("SplitLocusIntervals interval (start=%s) has size %d outside of bounds ideal=%d, max=%d",
                            split.get(0), splitSize, ideal, maxSize));
        }

        // test that every base in original is covered once by a base in split by locus intervals
        String diff = IntervalUtils.equateIntervals(test.originalIntervals, flat);
        Assert.assertNull(diff, diff);
    }

    /** small scale tests where the expected cuts are enumerated upfront for testing */
    private class SplitLocusIntervalsSmallTest extends TestDataProvider {
        final List<GenomeLoc> original;
        final public int parts;
        final public int expectedParts;
        final List<GenomeLoc> expected;

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected) {
            this(name, originalIntervals, parts,  expected, parts);
        }

        private SplitLocusIntervalsSmallTest(final String name, List<GenomeLoc> originalIntervals, final int parts, List<GenomeLoc> expected, int expectedParts) {
            super(SplitLocusIntervalsSmallTest.class, name);
            this.parts = parts;
            this.expectedParts = expectedParts;
            this.original = originalIntervals;
            this.expected = expected;
        }

        public String toString() {
            return String.format("%s parts=%d", super.toString(), parts);
        }
    }

    @DataProvider(name = "SplitLocusIntervalsSmallTest")
    public Object[][] createSplitLocusIntervalsSmallTest() {
        GenomeLoc bp01_10 = hg19GenomeLocParser.createGenomeLoc("1", 1, 10);

        GenomeLoc bp1_5 = hg19GenomeLocParser.createGenomeLoc("1", 1, 5);
        GenomeLoc bp6_10 = hg19GenomeLocParser.createGenomeLoc("1", 6, 10);
        new SplitLocusIntervalsSmallTest("cut into two", Arrays.asList(bp01_10), 2, Arrays.asList(bp1_5, bp6_10));

        GenomeLoc bp20_30 = hg19GenomeLocParser.createGenomeLoc("1", 20, 30);
        new SplitLocusIntervalsSmallTest("two in two", Arrays.asList(bp01_10, bp20_30), 2, Arrays.asList(bp01_10, bp20_30));

        GenomeLoc bp1_7 = hg19GenomeLocParser.createGenomeLoc("1", 1, 7);
        GenomeLoc bp8_10 = hg19GenomeLocParser.createGenomeLoc("1", 8, 10);
        GenomeLoc bp20_23 = hg19GenomeLocParser.createGenomeLoc("1", 20, 23);
        GenomeLoc bp24_30 = hg19GenomeLocParser.createGenomeLoc("1", 24, 30);
        new SplitLocusIntervalsSmallTest("two in three", Arrays.asList(bp01_10, bp20_30), 3,
                Arrays.asList(bp1_7, bp8_10, bp20_23, bp24_30));

        GenomeLoc bp1_2 = hg19GenomeLocParser.createGenomeLoc("1", 1, 2);
        GenomeLoc bp1_1 = hg19GenomeLocParser.createGenomeLoc("1", 1, 1);
        GenomeLoc bp2_2 = hg19GenomeLocParser.createGenomeLoc("1", 2, 2);
        new SplitLocusIntervalsSmallTest("too many pieces", Arrays.asList(bp1_2), 5, Arrays.asList(bp1_1, bp2_2), 2);

        new SplitLocusIntervalsSmallTest("emptyList", Collections.<GenomeLoc>emptyList(), 5, Collections.<GenomeLoc>emptyList(), 0);

        return SplitLocusIntervalsSmallTest.getTests(SplitLocusIntervalsSmallTest.class);
    }

    @Test(enabled = true, dataProvider = "SplitLocusIntervalsSmallTest")
    public void splitLocusIntervalsSmallTest(SplitLocusIntervalsSmallTest test) {
        List<List<GenomeLoc>> splitByLocus = IntervalUtils.splitLocusIntervals(test.original, test.parts);
        Assert.assertEquals(splitByLocus.size(), test.expectedParts, "SplitLocusIntervals failed to generate correct number of intervals");
        List<GenomeLoc> flat = IntervalUtils.flattenSplitIntervals(splitByLocus);

        // test sizes
        final long originalSize = IntervalUtils.intervalSize(test.original);
        final long splitSize = IntervalUtils.intervalSize(flat);
        Assert.assertEquals(splitSize, originalSize, "SplitLocusIntervals locs cover an incorrect number of bases");

        Assert.assertEquals(flat, test.expected, "SplitLocusIntervals locs not expected intervals");
    }

    //
    // Misc. tests
    //

    @Test(expectedExceptions=UserException.class)
    public void testMergeListsBySetOperatorNoOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> listEveryTwoFromOne = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, null);
        Assert.assertEquals(ret.size(), 100);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 0);
    }

    @Test
    public void testMergeListsBySetOperatorAllOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
            allSites.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 150);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 50);
    }

    @Test
    public void testMergeListsBySetOperator() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 5 == 0) {
                listEveryTwoFromTwo.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
                allSites.add(hg18GenomeLocParser.createGenomeLoc("chr1",x,x));
            }
        }

        List<GenomeLoc> ret;
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, null);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 20);
    }

    @Test
    public void testOverlappingIntervalsFromSameSourceWithIntersection() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> source1 = new ArrayList<GenomeLoc>();
        List<GenomeLoc> source2 = new ArrayList<GenomeLoc>();

        source1.add(hg18GenomeLocParser.createGenomeLoc("chr1", 10, 20));
        source1.add(hg18GenomeLocParser.createGenomeLoc("chr1", 15, 25));

        source2.add(hg18GenomeLocParser.createGenomeLoc("chr1", 16, 18));
        source2.add(hg18GenomeLocParser.createGenomeLoc("chr1", 22, 24));

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(source1, source2, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 2);
    }

    @Test
    public void testGetContigLengths() {
        Map<String, Integer> lengths = IntervalUtils.getContigSizes(new File(BaseTest.hg18Reference));
        Assert.assertEquals((long)lengths.get("chr1"), 247249719);
        Assert.assertEquals((long)lengths.get("chr2"), 242951149);
        Assert.assertEquals((long)lengths.get("chr3"), 199501827);
        Assert.assertEquals((long)lengths.get("chr20"), 62435964);
        Assert.assertEquals((long)lengths.get("chrX"), 154913754);
    }

    @Test
    public void testParseIntervalArguments() {
        Assert.assertEquals(getLocs().size(), 45);
        Assert.assertEquals(getLocs("chr1", "chr2", "chr3").size(), 3);
        Assert.assertEquals(getLocs("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2").size(), 4);
    }

    @Test
    public void testIsIntervalFile() {
        Assert.assertTrue(IntervalUtils.isIntervalFile(BaseTest.privateTestDir + "empty_intervals.list"));
        Assert.assertTrue(IntervalUtils.isIntervalFile(BaseTest.privateTestDir + "empty_intervals.list", true));

        List<String> extensions = Arrays.asList("bed", "interval_list", "intervals", "list", "picard");
        for (String extension: extensions) {
            Assert.assertTrue(IntervalUtils.isIntervalFile("test_intervals." + extension, false), "Tested interval file extension: " + extension);
        }
    }

    @Test(expectedExceptions = UserException.CouldNotReadInputFile.class)
    public void testMissingIntervalFile() {
        IntervalUtils.isIntervalFile(BaseTest.privateTestDir + "no_such_intervals.list");
    }

    @Test
    public void testFixedScatterIntervalsBasic() {
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3");

        List<File> files = testFiles("basic.", 3, ".intervals");

        List<GenomeLoc> locs = getLocs("chr1", "chr2", "chr3");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterFixedIntervalsLessFiles() {
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3");
        GenomeLoc chr4 = hg18GenomeLocParser.parseGenomeLoc("chr4");

        List<File> files = testFiles("less.", 3, ".intervals");

        List<GenomeLoc> locs = getLocs("chr1", "chr2", "chr3", "chr4");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
        Assert.assertEquals(locs3.get(1), chr4);
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testSplitFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = getLocs("chr1", "chr2");
        IntervalUtils.splitFixedIntervals(locs, files.size());
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testScatterFixedIntervalsMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        List<GenomeLoc> locs = getLocs("chr1", "chr2");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, locs.size()); // locs.size() instead of files.size()
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);
    }
    @Test
    public void testScatterFixedIntervalsStart() {
        List<String> intervals = Arrays.asList("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2");
        GenomeLoc chr1a = hg18GenomeLocParser.parseGenomeLoc("chr1:1-2");
        GenomeLoc chr1b = hg18GenomeLocParser.parseGenomeLoc("chr1:4-5");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:1-1");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = getLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs2.get(0), chr1b);
        Assert.assertEquals(locs3.get(0), chr2);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsMiddle() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:1-2", "chr2:4-5", "chr3:2-2");
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1");
        GenomeLoc chr2a = hg18GenomeLocParser.parseGenomeLoc("chr2:1-2");
        GenomeLoc chr2b = hg18GenomeLocParser.parseGenomeLoc("chr2:4-5");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = getLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs3.get(0), chr2b);
        Assert.assertEquals(locs3.get(1), chr3);
    }

    @Test
    public void testScatterFixedIntervalsEnd() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:2-2", "chr3:1-2", "chr3:4-5");
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:2-2");
        GenomeLoc chr3a = hg18GenomeLocParser.parseGenomeLoc("chr3:1-2");
        GenomeLoc chr3b = hg18GenomeLocParser.parseGenomeLoc("chr3:4-5");

        List<File> files = testFiles("split.", 3, ".intervals");

        List<GenomeLoc> locs = getLocs(intervals);
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());
        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3a);
        Assert.assertEquals(locs3.get(0), chr3b);
    }

    @Test
    public void testScatterFixedIntervalsFile() {
        List<File> files = testFiles("sg.", 20, ".intervals");
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(BaseTest.GATKDataLocation + "whole_exome_agilent_designed_120.targets.hg18.chr20.interval_list"));
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(locs, files.size());

        int[] counts = {
                125, 138, 287, 291, 312, 105, 155, 324,
                295, 298, 141, 121, 285, 302, 282, 88,
                116, 274, 282, 248
//                5169, 5573, 10017, 10567, 10551,
//                5087, 4908, 10120, 10435, 10399,
//                5391, 4735, 10621, 10352, 10654,
//                5227, 5256, 10151, 9649, 9825
        };

        //String splitCounts = "";
        for (int i = 0; i < splits.size(); i++) {
            int splitCount = splits.get(i).size();
            Assert.assertEquals(splitCount, counts[i], "Num intervals in split " + i);
        }
        //System.out.println(splitCounts.substring(2));

        IntervalUtils.scatterFixedIntervals(hg18Header, splits, files);

        int locIndex = 0;
        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), counts[i], "Intervals in " + file);
            for (GenomeLoc parsedLoc: parsedLocs)
                Assert.assertEquals(parsedLoc, locs.get(locIndex), String.format("Genome loc %d from file %d", locIndex++, i));
        }
        Assert.assertEquals(locIndex, locs.size(), "Total number of GenomeLocs");
    }

    @Test
    public void testScatterFixedIntervalsMax() {
        List<File> files = testFiles("sg.", 85, ".intervals");
        List<List<GenomeLoc>> splits = IntervalUtils.splitFixedIntervals(hg19ReferenceLocs, files.size());
        IntervalUtils.scatterFixedIntervals(hg19Header, splits, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    @Test
    public void testScatterContigIntervalsOrder() {
        List<String> intervals = Arrays.asList("chr2:1-1", "chr1:1-1", "chr3:2-2");
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:1-1");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr2);
        Assert.assertEquals(locs2.get(0), chr1);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsBasic() {
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3");

        List<File> files = testFiles("contig_basic.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs("chr1", "chr2", "chr3"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsLessFiles() {
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3");
        GenomeLoc chr4 = hg18GenomeLocParser.parseGenomeLoc("chr4");

        List<File> files = testFiles("contig_less.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs("chr1", "chr2", "chr3", "chr4"), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3);
        Assert.assertEquals(locs3.get(0), chr4);
    }

    @Test(expectedExceptions=UserException.BadInput.class)
    public void testScatterContigIntervalsMoreFiles() {
        List<File> files = testFiles("contig_more.", 3, ".intervals");
        IntervalUtils.scatterContigIntervals(hg18Header, getLocs("chr1", "chr2"), files);
    }

    @Test
    public void testScatterContigIntervalsStart() {
        List<String> intervals = Arrays.asList("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2");
        GenomeLoc chr1a = hg18GenomeLocParser.parseGenomeLoc("chr1:1-2");
        GenomeLoc chr1b = hg18GenomeLocParser.parseGenomeLoc("chr1:4-5");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:1-1");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:2-2");

        List<File> files = testFiles("contig_split_start.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs1.get(1), chr1b);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsMiddle() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:1-2", "chr2:4-5", "chr3:2-2");
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1");
        GenomeLoc chr2a = hg18GenomeLocParser.parseGenomeLoc("chr2:1-2");
        GenomeLoc chr2b = hg18GenomeLocParser.parseGenomeLoc("chr2:4-5");
        GenomeLoc chr3 = hg18GenomeLocParser.parseGenomeLoc("chr3:2-2");

        List<File> files = testFiles("contig_split_middle.", 3, ".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 2);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs2.get(1), chr2b);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterContigIntervalsEnd() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:2-2", "chr3:1-2", "chr3:4-5");
        GenomeLoc chr1 = hg18GenomeLocParser.parseGenomeLoc("chr1:1-1");
        GenomeLoc chr2 = hg18GenomeLocParser.parseGenomeLoc("chr2:2-2");
        GenomeLoc chr3a = hg18GenomeLocParser.parseGenomeLoc("chr3:1-2");
        GenomeLoc chr3b = hg18GenomeLocParser.parseGenomeLoc("chr3:4-5");

        List<File> files = testFiles("contig_split_end.", 3 ,".intervals");

        IntervalUtils.scatterContigIntervals(hg18Header, getLocs(intervals), files);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(0).toString()));
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(1).toString()));
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Arrays.asList(files.get(2).toString()));

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3a);
        Assert.assertEquals(locs3.get(1), chr3b);
    }

    @Test
    public void testScatterContigIntervalsMax() {
        List<File> files = testFiles("sg.", 85, ".intervals");
        IntervalUtils.scatterContigIntervals(hg19Header, hg19ReferenceLocs, files);

        for (int i = 0; i < files.size(); i++) {
            String file = files.get(i).toString();
            List<GenomeLoc> parsedLocs = IntervalUtils.parseIntervalArguments(hg19GenomeLocParser, Arrays.asList(file));
            Assert.assertEquals(parsedLocs.size(), 1, "parsedLocs[" + i + "].size()");
            Assert.assertEquals(parsedLocs.get(0), hg19ReferenceLocs.get(i), "parsedLocs[" + i + "].get()");
        }
    }

    private List<File> testFiles(String prefix, int count, String suffix) {
        ArrayList<File> files = new ArrayList<File>();
        for (int i = 1; i <= count; i++) {
            files.add(createTempFile(prefix + i, suffix));
        }
        return files;
    }

    @DataProvider(name="unmergedIntervals")
    public Object[][] getUnmergedIntervals() {
        return new Object[][] {
                new Object[] {"small_unmerged_picard_intervals.list"},
                new Object[] {"small_unmerged_gatk_intervals.list"}
        };
    }

    @Test(dataProvider="unmergedIntervals")
    public void testUnmergedIntervals(String unmergedIntervals) {
        List<GenomeLoc> locs = IntervalUtils.parseIntervalArguments(hg18GenomeLocParser, Collections.singletonList(privateTestDir + unmergedIntervals));
        Assert.assertEquals(locs.size(), 2);

        List<GenomeLoc> merged;

        merged = IntervalUtils.mergeIntervalLocations(locs, IntervalMergingRule.ALL);
        Assert.assertEquals(merged.size(), 1);

        // Test that null means the same as ALL
        merged = IntervalUtils.mergeIntervalLocations(locs, null);
        Assert.assertEquals(merged.size(), 1);
    }

    /*
    Split into tests that can be written to files and tested by writeFlankingIntervals,
    and lists that cannot but are still handled by getFlankingIntervals.
    */
    private static abstract class FlankingIntervalsTestData extends TestDataProvider {
        final public File referenceFile;
        final public GenomeLocParser parser;
        final int basePairs;
        final List<GenomeLoc> original;
        final List<GenomeLoc> expected;

        protected FlankingIntervalsTestData(Class<?> clazz, String name, File referenceFile, GenomeLocParser parser,
                                          int basePairs, List<String> original, List<String> expected) {
            super(clazz, name);
            this.referenceFile = referenceFile;
            this.parser = parser;
            this.basePairs = basePairs;
            this.original = parse(parser, original);
            this.expected = parse(parser, expected);
        }

        private static List<GenomeLoc> parse(GenomeLocParser parser, List<String> locs) {
            List<GenomeLoc> parsed = new ArrayList<GenomeLoc>();
            for (String loc: locs)
                parsed.add("unmapped".equals(loc) ? GenomeLoc.UNMAPPED : parser.parseGenomeLoc(loc));
            return parsed;
        }
    }

    private static class FlankingIntervalsFile extends FlankingIntervalsTestData {
        public FlankingIntervalsFile(String name, File referenceFile, GenomeLocParser parser,
                                     int basePairs, List<String> original, List<String> expected) {
            super(FlankingIntervalsFile.class, name, referenceFile, parser, basePairs, original, expected);
        }
    }

    private static class FlankingIntervalsList extends FlankingIntervalsTestData {
        public FlankingIntervalsList(String name, File referenceFile, GenomeLocParser parser,
                                     int basePairs, List<String> original, List<String> expected) {
            super(FlankingIntervalsList.class, name, referenceFile, parser, basePairs, original, expected);
        }
    }

    /* Intervals where the original and the flanks can be written to files. */
    @DataProvider(name = "flankingIntervalsFiles")
    public Object[][] getFlankingIntervalsFiles() {
        File hg19ReferenceFile = new File(BaseTest.hg19Reference);
        int hg19Length1 = hg19GenomeLocParser.getContigInfo("1").getSequenceLength();

        new FlankingIntervalsFile("atStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:1"),
                Arrays.asList("1:2"));

        new FlankingIntervalsFile("atStartBase50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:1"),
                Arrays.asList("1:2-51"));

        new FlankingIntervalsFile("atStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:1-10"),
                Arrays.asList("1:11-60"));

        new FlankingIntervalsFile("atEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + hg19Length1),
                Arrays.asList("1:" + (hg19Length1 - 1)));

        new FlankingIntervalsFile("atEndBase50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:" + hg19Length1),
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 50, hg19Length1 - 1)));

        new FlankingIntervalsFile("atEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 10, hg19Length1)),
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 60, hg19Length1 - 11)));

        new FlankingIntervalsFile("nearStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:2"),
                Arrays.asList("1:1", "1:3"));

        new FlankingIntervalsFile("nearStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:21-30"),
                Arrays.asList("1:1-20", "1:31-80"));

        new FlankingIntervalsFile("nearEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + (hg19Length1 - 1)),
                Arrays.asList("1:" + (hg19Length1 - 2), "1:" + hg19Length1));

        new FlankingIntervalsFile("nearEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 30, hg19Length1 - 21)),
                Arrays.asList(
                        String.format("1:%d-%d", hg19Length1 - 80, hg19Length1 - 31),
                        String.format("1:%d-%d", hg19Length1 - 20, hg19Length1)));

        new FlankingIntervalsFile("beyondStartBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:3"),
                Arrays.asList("1:2", "1:4"));

        new FlankingIntervalsFile("beyondStartRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200"),
                Arrays.asList("1:51-100", "1:201-250"));

        new FlankingIntervalsFile("beyondEndBase1", hg19ReferenceFile, hg19GenomeLocParser, 1,
                Arrays.asList("1:" + (hg19Length1 - 3)),
                Arrays.asList("1:" + (hg19Length1 - 4), "1:" + (hg19Length1 - 2)));

        new FlankingIntervalsFile("beyondEndRange50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList(String.format("1:%d-%d", hg19Length1 - 200, hg19Length1 - 101)),
                Arrays.asList(
                        String.format("1:%d-%d", hg19Length1 - 250, hg19Length1 - 201),
                        String.format("1:%d-%d", hg19Length1 - 100, hg19Length1 - 51)));

        new FlankingIntervalsFile("betweenFar50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:401-500"),
                Arrays.asList("1:51-100", "1:201-250", "1:351-400", "1:501-550"));

        new FlankingIntervalsFile("betweenSpan50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:301-400"),
                Arrays.asList("1:51-100", "1:201-300", "1:401-450"));

        new FlankingIntervalsFile("betweenOverlap50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:271-400"),
                Arrays.asList("1:51-100", "1:201-270", "1:401-450"));

        new FlankingIntervalsFile("betweenShort50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:221-400"),
                Arrays.asList("1:51-100", "1:201-220", "1:401-450"));

        new FlankingIntervalsFile("betweenNone50", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:121-400"),
                Arrays.asList("1:51-100", "1:401-450"));

        new FlankingIntervalsFile("twoContigs", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "2:301-400"),
                Arrays.asList("1:51-100", "1:201-250", "2:251-300", "2:401-450"));

        // Explicit testing a problematic agilent target pair
        new FlankingIntervalsFile("badAgilent", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("2:74756257-74756411", "2:74756487-74756628"),
                // wrong!    ("2:74756206-74756256", "2:74756412-74756462", "2:74756436-74756486", "2:74756629-74756679")
                Arrays.asList("2:74756207-74756256", "2:74756412-74756486", "2:74756629-74756678"));

        return TestDataProvider.getTests(FlankingIntervalsFile.class);
    }

    /* Intervals where either the original and/or the flanks cannot be written to a file. */
    @DataProvider(name = "flankingIntervalsLists")
    public Object[][] getFlankingIntervalsLists() {
        File hg19ReferenceFile = new File(BaseTest.hg19Reference);
        List<String> empty = Collections.emptyList();

        new FlankingIntervalsList("empty", hg19ReferenceFile, hg19GenomeLocParser, 50,
                empty,
                empty);

        new FlankingIntervalsList("unmapped", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("unmapped"),
                empty);

        new FlankingIntervalsList("fullContig", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1"),
                empty);

        new FlankingIntervalsList("fullContigs", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1", "2", "3"),
                empty);

        new FlankingIntervalsList("betweenWithUnmapped", hg19ReferenceFile, hg19GenomeLocParser, 50,
                Arrays.asList("1:101-200", "1:301-400", "unmapped"),
                Arrays.asList("1:51-100", "1:201-300", "1:401-450"));

        return TestDataProvider.getTests(FlankingIntervalsList.class);
    }

    @Test(dataProvider = "flankingIntervalsFiles")
    public void testWriteFlankingIntervals(FlankingIntervalsTestData data) throws Exception {
        File originalFile = createTempFile("original.", ".intervals");
        File flankingFile = createTempFile("flanking.", ".intervals");
        try {
            List<String> lines = new ArrayList<String>();
            for (GenomeLoc loc: data.original)
                lines.add(loc.toString());
            FileUtils.writeLines(originalFile, lines);

            IntervalUtils.writeFlankingIntervals(data.referenceFile, originalFile, flankingFile, data.basePairs);

            List<GenomeLoc> actual = IntervalUtils.intervalFileToList(data.parser, flankingFile.getAbsolutePath());

            String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                    data.toString(), data.original, actual, data.expected);
            Assert.assertEquals(actual, data.expected, description);
        } finally {
            FileUtils.deleteQuietly(originalFile);
            FileUtils.deleteQuietly(flankingFile);
        }
    }

    @Test(dataProvider = "flankingIntervalsLists", expectedExceptions = UserException.class)
    public void testWritingBadFlankingIntervals(FlankingIntervalsTestData data) throws Exception {
        File originalFile = createTempFile("original.", ".intervals");
        File flankingFile = createTempFile("flanking.", ".intervals");
        try {
            List<String> lines = new ArrayList<String>();
            for (GenomeLoc loc: data.original)
                lines.add(loc.toString());
            FileUtils.writeLines(originalFile, lines);

            // Should throw a user exception on bad input if either the original
            // intervals are empty or if the flanking intervals are empty
            IntervalUtils.writeFlankingIntervals(data.referenceFile, originalFile, flankingFile, data.basePairs);
        } finally {
            FileUtils.deleteQuietly(originalFile);
            FileUtils.deleteQuietly(flankingFile);
        }
    }

    @Test(dataProvider = "flankingIntervalsLists")
    public void testGetFlankingIntervals(FlankingIntervalsTestData data) {
        List<GenomeLoc> actual = IntervalUtils.getFlankingIntervals(data.parser, data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.expected);
        Assert.assertEquals(actual, data.expected, description);
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testExceptionUponLegacyIntervalSyntax() throws Exception {
        final GenomeLocParser parser = new GenomeLocParser(new CachingIndexedFastaSequenceFile(new File(BaseTest.hg19Reference)));

        // Attempting to use the legacy -L "interval1;interval2" syntax should produce an exception:
        IntervalBinding<Feature> binding = new IntervalBinding<Feature>("1;2");
        binding.getIntervals(parser);
    }

    @DataProvider(name="invalidIntervalTestData")
    public Object[][] invalidIntervalDataProvider() throws Exception {
        File fastaFile = new File(publicTestDir + "exampleFASTA.fasta");
        GenomeLocParser genomeLocParser = new GenomeLocParser(new IndexedFastaSequenceFile(fastaFile));

        return new Object[][] {
                new Object[] {genomeLocParser, "chr1", 10000000, 20000000},
                new Object[] {genomeLocParser, "chr2", 1, 2},
                new Object[] {genomeLocParser, "chr1", -1, 50}
        };
    }

    /*
     * This test is disabled because its assumption that we will not throw an error
     * upon parsing invalid Picard intervals is no longer true, as htsjdk has added
     * extra protection against invalid intervals to IntervalList.add().
     *
     * We should reconsider our decision in IntervalUtils.intervalFileToList() to
     * silently ignore invalid intervals when parsing Picard interval files, as it's
     * inconsistent with the way we handle invalid intervals for GATK interval files
     * (throw a UserException, covered by testInvalidGATKFileIntervalHandling() below),
     * and update this test accordingly.
     */
    @Test(dataProvider="invalidIntervalTestData", enabled = false)
    public void testInvalidPicardIntervalHandling(GenomeLocParser genomeLocParser,
                                                  String contig, int intervalStart, int intervalEnd ) throws Exception {

        SAMFileHeader picardFileHeader = new SAMFileHeader();
        picardFileHeader.addSequence(genomeLocParser.getContigInfo("chr1"));
        IntervalList picardIntervals = new IntervalList(picardFileHeader);
        picardIntervals.add(new Interval(contig, intervalStart, intervalEnd, true, "dummyname"));

        File picardIntervalFile = createTempFile("testInvalidPicardIntervalHandling", ".intervals");
        picardIntervals.write(picardIntervalFile);

        List<IntervalBinding<Feature>> intervalArgs = new ArrayList<IntervalBinding<Feature>>(1);
        intervalArgs.add(new IntervalBinding<Feature>(picardIntervalFile.getAbsolutePath()));

        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    @Test(expectedExceptions=UserException.class, dataProvider="invalidIntervalTestData")
    public void testInvalidGATKFileIntervalHandling(GenomeLocParser genomeLocParser,
                                                    String contig, int intervalStart, int intervalEnd ) throws Exception {

        File gatkIntervalFile = createTempFile("testInvalidGATKFileIntervalHandling", ".intervals",
                String.format("%s:%d-%d", contig, intervalStart, intervalEnd));

        List<IntervalBinding<Feature>> intervalArgs = new ArrayList<IntervalBinding<Feature>>(1);
        intervalArgs.add(new IntervalBinding<Feature>(gatkIntervalFile.getAbsolutePath()));

        IntervalUtils.loadIntervals(intervalArgs, IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, genomeLocParser);
    }

    private File createTempFile( String tempFilePrefix, String tempFileExtension, String... lines ) throws Exception {
        File tempFile = BaseTest.createTempFile(tempFilePrefix, tempFileExtension);
        FileUtils.writeLines(tempFile, Arrays.asList(lines));
        return tempFile;
    }

    @DataProvider(name = "sortAndMergeIntervals")
    public Object[][] getSortAndMergeIntervals() {
        return new Object[][] {
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, getLocs("chr1:1", "chr1:3", "chr1:2"), getLocs("chr1:1", "chr1:2", "chr1:3") },
                new Object[] { IntervalMergingRule.ALL, getLocs("chr1:1", "chr1:3", "chr1:2"), getLocs("chr1:1-3") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, getLocs("chr1:1", "chr1:3", "chr2:2"), getLocs("chr1:1", "chr1:3", "chr2:2") },
                new Object[] { IntervalMergingRule.ALL, getLocs("chr1:1", "chr1:3", "chr2:2"), getLocs("chr1:1", "chr1:3", "chr2:2") },
                new Object[] { IntervalMergingRule.OVERLAPPING_ONLY, getLocs("chr1:1", "chr1"), getLocs("chr1") },
                new Object[] { IntervalMergingRule.ALL, getLocs("chr1:1", "chr1"), getLocs("chr1") }
        };
    }

    @Test(dataProvider = "sortAndMergeIntervals")
    public void testSortAndMergeIntervals(IntervalMergingRule merge, List<GenomeLoc> unsorted, List<GenomeLoc> expected) {
        List<GenomeLoc> sorted = IntervalUtils.sortAndMergeIntervals(hg18GenomeLocParser, unsorted, merge).toList();
        Assert.assertEquals(sorted, expected);
    }
}
