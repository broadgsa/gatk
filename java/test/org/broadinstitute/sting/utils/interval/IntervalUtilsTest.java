package org.broadinstitute.sting.utils.interval;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

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



    @BeforeClass
    public static void init() throws FileNotFoundException {
        seq = new IndexedFastaSequenceFile(new File(hg18Reference));
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    @Test
    public void testMergeListsBySetOperatorNoOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> listEveryTwoFromOne = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
            else
                listEveryTwoFromOne.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.UNION);
        Assert.assertEquals(100,ret.size());
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, listEveryTwoFromOne, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(0,ret.size());
    }

    @Test
    public void testMergeListsBySetOperatorAllOverlap() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 2 == 0)
                listEveryTwoFromTwo.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
            allSites.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(150,ret.size());
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(50,ret.size());
    }

    @Test
    public void testMergeListsBySetOperator() {
        // a couple of lists we'll use for the testing
        List<GenomeLoc> allSites = new ArrayList<GenomeLoc>();
        List<GenomeLoc> listEveryTwoFromTwo = new ArrayList<GenomeLoc>();

        // create the two lists we'll use
        for (int x = 1; x < 101; x++) {
            if (x % 5 == 0) {
                listEveryTwoFromTwo.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
                allSites.add(GenomeLocParser.createGenomeLoc("chr1",x,x));
            }
        }

        List<GenomeLoc> ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.UNION);
        Assert.assertEquals(40,ret.size());
        ret = IntervalUtils.mergeListsBySetOperator(listEveryTwoFromTwo, allSites, IntervalSetRule.INTERSECTION);
        Assert.assertEquals(20,ret.size());
    }
}
