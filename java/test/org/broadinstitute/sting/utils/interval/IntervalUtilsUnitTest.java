package org.broadinstitute.sting.utils.interval;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * test out the interval utility methods
 */
public class IntervalUtilsUnitTest extends BaseTest {
    // used to seed the genome loc parser with a sequence dictionary
    private static File reference = new File(BaseTest.hg18Reference);
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() {
        ReferenceSequenceFile seq = new IndexedFastaSequenceFile(reference);
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test
    public void testMergeListsBySetOperatorNoOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> listEveryTwoFromOne = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(genomeLocParser.createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(genomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.UNION);
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
                listEveryTwoFromTwo.add(genomeLocParser.createGenomeLoc("chr1",x,x));
            allSites.add(genomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
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
                listEveryTwoFromTwo.add(genomeLocParser.createGenomeLoc("chr1",x,x));
                allSites.add(genomeLocParser.createGenomeLoc("chr1",x,x));
            }
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(ret.size(), 40);
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(ret.size(), 20);
    }

    @Test
    public void testCountContigs() {
        List<String> chrs = new ArrayList<String>();
        for (int i = 1; i <= 22; i++)
            chrs.add("chr" + i);
        chrs.add("chrX");
        chrs.add("chrY");

        List<String> chrsNoRandom = Arrays.asList("chr12", "chr14", "chr20", "chrY");
        List<String> chrsWithRandom = new ArrayList<String>();
        chrsWithRandom.add("chrM");
        chrsWithRandom.addAll(chrs);
        for (String chr: chrs)
            if(!chrsNoRandom.contains(chr))
                chrsWithRandom.add(chr + "_random");

        Assert.assertEquals(IntervalUtils.distinctContigs(reference), chrsWithRandom);
        Assert.assertEquals(IntervalUtils.distinctContigs(reference, Arrays.asList(BaseTest.validationDataLocation + "TCGA-06-0188.interval_list")), chrs);
        Assert.assertEquals(IntervalUtils.distinctContigs(reference, Arrays.asList("chr1:1-1", "chr2:1-1", "chr3:2-2")), Arrays.asList("chr1","chr2","chr3"));
        Assert.assertEquals(IntervalUtils.distinctContigs(reference, Arrays.asList("chr2:1-1", "chr1:1-1", "chr3:2-2")), Arrays.asList("chr1","chr2","chr3"));
    }

    @Test
    public void testBasicScatter() {
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3");

        List<File> files = testFiles("basic.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2", "chr3"), files, false);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterLessFiles() {
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3");
        GenomeLoc chr4 = genomeLocParser.parseGenomeInterval("chr4");

        List<File> files = testFiles("less.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2", "chr3", "chr4"), files, false);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs1.get(1), chr2);
        Assert.assertEquals(locs2.get(0), chr3);
        Assert.assertEquals(locs3.get(0), chr4);
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testScatterMoreFiles() {
        List<File> files = testFiles("more.", 3, ".intervals");
        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2"), files, false);
    }

    @Test
    public void testScatterIntervals() {
        List<String> intervals = Arrays.asList("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2");
        GenomeLoc chr1a = genomeLocParser.parseGenomeInterval("chr1:1-2");
        GenomeLoc chr1b = genomeLocParser.parseGenomeInterval("chr1:4-5");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2:1-1");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, intervals, files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs1.get(1), chr1b);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test(enabled=false) // disabled, GenomeLoc.compareTo() returns 0 for two locs with the same start, causing an exception in GLSS.add().
    public void testScatterIntervalsWithTheSameStart() {
        List<File> files = testFiles("sg.", 20, ".intervals");
        IntervalUtils.scatterIntervalArguments(new File(hg18Reference), Arrays.asList(BaseTest.GATKDataLocation + "whole_exome_agilent_designed_120.targets.hg18.chr20.interval_list"), files, false);
    }

    @Test
    public void testScatterOrder() {
        List<String> intervals = Arrays.asList("chr2:1-1", "chr1:1-1", "chr3:2-2");
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1:1-1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2:1-1");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3:2-2");

        List<File> files = testFiles("split.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, intervals, files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testBasicScatterByContig() {
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3");

        List<File> files = testFiles("contig_basic.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2", "chr3"), files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterByContigLessFiles() {
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3");
        GenomeLoc chr4 = genomeLocParser.parseGenomeInterval("chr4");

        List<File> files = testFiles("contig_less.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2", "chr3", "chr4"), files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
        Assert.assertEquals(locs3.get(1), chr4);
    }

    @Test(expectedExceptions=UserException.BadArgumentValue.class)
    public void testScatterByContigMoreFiles() {
        List<File> files = testFiles("contig_more.", 3, ".intervals");
        IntervalUtils.scatterIntervalArguments(reference, Arrays.asList("chr1", "chr2"), files, true);
    }

    @Test
    public void testScatterByContigIntervalsStart() {
        List<String> intervals = Arrays.asList("chr1:1-2", "chr1:4-5", "chr2:1-1", "chr3:2-2");
        GenomeLoc chr1a = genomeLocParser.parseGenomeInterval("chr1:1-2");
        GenomeLoc chr1b = genomeLocParser.parseGenomeInterval("chr1:4-5");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2:1-1");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3:2-2");

        List<File> files = testFiles("contig_split_start.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, intervals, files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 2);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1a);
        Assert.assertEquals(locs1.get(1), chr1b);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterByContigIntervalsMiddle() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:1-2", "chr2:4-5", "chr3:2-2");
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1:1-1");
        GenomeLoc chr2a = genomeLocParser.parseGenomeInterval("chr2:1-2");
        GenomeLoc chr2b = genomeLocParser.parseGenomeInterval("chr2:4-5");
        GenomeLoc chr3 = genomeLocParser.parseGenomeInterval("chr3:2-2");

        List<File> files = testFiles("contig_split_middle.", 3, ".intervals");

        IntervalUtils.scatterIntervalArguments(reference, intervals, files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 2);
        Assert.assertEquals(locs3.size(), 1);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2a);
        Assert.assertEquals(locs2.get(1), chr2b);
        Assert.assertEquals(locs3.get(0), chr3);
    }

    @Test
    public void testScatterByContigIntervalsEnd() {
        List<String> intervals = Arrays.asList("chr1:1-1", "chr2:2-2", "chr3:1-2", "chr3:4-5");
        GenomeLoc chr1 = genomeLocParser.parseGenomeInterval("chr1:1-1");
        GenomeLoc chr2 = genomeLocParser.parseGenomeInterval("chr2:2-2");
        GenomeLoc chr3a = genomeLocParser.parseGenomeInterval("chr3:1-2");
        GenomeLoc chr3b = genomeLocParser.parseGenomeInterval("chr3:4-5");

        List<File> files = testFiles("contig_split_end.", 3 ,".intervals");

        IntervalUtils.scatterIntervalArguments(reference, intervals, files, true);

        List<GenomeLoc> locs1 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(0).toString()), false);
        List<GenomeLoc> locs2 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(1).toString()), false);
        List<GenomeLoc> locs3 = IntervalUtils.parseIntervalArguments(genomeLocParser, Arrays.asList(files.get(2).toString()), false);

        Assert.assertEquals(locs1.size(), 1);
        Assert.assertEquals(locs2.size(), 1);
        Assert.assertEquals(locs3.size(), 2);

        Assert.assertEquals(locs1.get(0), chr1);
        Assert.assertEquals(locs2.get(0), chr2);
        Assert.assertEquals(locs3.get(0), chr3a);
        Assert.assertEquals(locs3.get(1), chr3b);
    }

    private List<File> testFiles(String prefix, int count, String suffix) {
        try {
            ArrayList<File> files = new ArrayList<File>();
            for (int i = 1; i <= count; i++) {
                File tmpFile = File.createTempFile(prefix + i, suffix);
                tmpFile.deleteOnExit();
                files.add(tmpFile);
            }
            return files;
        } catch (IOException e) {
            throw new UserException.BadTmpDir("Unable to create temp file: " + e);
        }
    }
}
