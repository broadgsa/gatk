// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;

import java.io.File;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocTest extends BaseTest {
    private static FastaSequenceFile2 seq;

    @BeforeClass
    public static void init() {
        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testIsBetween() {
        logger.warn("Executing testIsBetween");
        GenomeLoc.setupRefContigOrdering(seq);
        GenomeLoc locMiddle = new GenomeLoc("chr1", 3, 3);

        GenomeLoc locLeft = new GenomeLoc("chr1", 1, 1);
        GenomeLoc locRight = new GenomeLoc("chr1", 5, 5);

        Assert.assertTrue(locMiddle.isBetween(locLeft, locRight));
        Assert.assertFalse(locLeft.isBetween(locMiddle, locRight));
        Assert.assertFalse(locRight.isBetween(locLeft, locMiddle));

    }
    @Test
    public void testContigIndex() {
        logger.warn("Executing testContigIndex");
        GenomeLoc locOne = new GenomeLoc("chr1",1,1);
        Assert.assertEquals(locOne.getContigIndex(), 1);
        Assert.assertEquals(locOne.getContig(), "chr1");

        GenomeLoc locX = new GenomeLoc("chrX",1,1);
        Assert.assertEquals(locX.getContigIndex(), 23);
        Assert.assertEquals(locX.getContig(), "chrX");

        GenomeLoc locNumber = new GenomeLoc(1,1,1);
        Assert.assertEquals(locNumber.getContigIndex(), 1);
        Assert.assertEquals(locNumber.getContig(), "chr1");
        Assert.assertEquals(locOne.compareTo(locNumber), 0);

    }

    @Test
    public void testCompareTo() {
        logger.warn("Executing testCompareTo");
        GenomeLoc twoOne = new GenomeLoc("chr2", 1);
        GenomeLoc twoFive = new GenomeLoc("chr2", 5);
        GenomeLoc twoOtherFive = new GenomeLoc("chr2", 5);
        Assert.assertEquals(0, twoFive.compareTo(twoOtherFive));

        Assert.assertEquals(-1, twoOne.compareTo(twoFive));
        Assert.assertEquals(1, twoFive.compareTo(twoOne));

        GenomeLoc oneOne = new GenomeLoc("chr1", 5);
        Assert.assertEquals(-1, oneOne.compareTo(twoOne));
        Assert.assertEquals(1, twoOne.compareTo(oneOne));
    }



}