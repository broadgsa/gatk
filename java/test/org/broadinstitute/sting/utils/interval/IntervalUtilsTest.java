package org.broadinstitute.sting.utils.interval;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * test out the interval utility methods
 */
public class IntervalUtilsTest extends BaseTest {
    // used to seed the genome loc parser with a sequence dictionary
    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        seq = new IndexedFastaSequenceFile(new File(hg18Reference));
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
}
