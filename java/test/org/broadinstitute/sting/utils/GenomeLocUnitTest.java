// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.broadinstitute.sting.BaseTest;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocUnitTest extends BaseTest {
    private static ReferenceSequenceFile seq;

    @BeforeClass
    public static void init() throws FileNotFoundException {
        // sequence
        seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testIsBetween() {
        logger.warn("Executing testIsBetween");

        GenomeLoc locMiddle = GenomeLocParser.createGenomeLoc("chr1", 3, 3);

        GenomeLoc locLeft = GenomeLocParser.createGenomeLoc("chr1", 1, 1);
        GenomeLoc locRight = GenomeLocParser.createGenomeLoc("chr1", 5, 5);

        Assert.assertTrue(locMiddle.isBetween(locLeft, locRight));
        Assert.assertFalse(locLeft.isBetween(locMiddle, locRight));
        Assert.assertFalse(locRight.isBetween(locLeft, locMiddle));

    }
    @Test
    public void testContigIndex() {
        logger.warn("Executing testContigIndex");
        GenomeLoc locOne = GenomeLocParser.createGenomeLoc("chr1",1,1);
        Assert.assertEquals(locOne.getContigIndex(), 1);
        Assert.assertEquals(locOne.getContig(), "chr1");

        GenomeLoc locX = GenomeLocParser.createGenomeLoc("chrX",1,1);
        Assert.assertEquals(locX.getContigIndex(), 23);
        Assert.assertEquals(locX.getContig(), "chrX");

        GenomeLoc locNumber = GenomeLocParser.createGenomeLoc(1,1,1);
        Assert.assertEquals(locNumber.getContigIndex(), 1);
        Assert.assertEquals(locNumber.getContig(), "chr1");
        Assert.assertEquals(locOne.compareTo(locNumber), 0);

    }

    @Test
    public void testCompareTo() {
        logger.warn("Executing testCompareTo");
        GenomeLoc twoOne = GenomeLocParser.createGenomeLoc("chr2", 1);
        GenomeLoc twoFive = GenomeLocParser.createGenomeLoc("chr2", 5);
        GenomeLoc twoOtherFive = GenomeLocParser.createGenomeLoc("chr2", 5);
        Assert.assertEquals(0, twoFive.compareTo(twoOtherFive));

        Assert.assertEquals(-1, twoOne.compareTo(twoFive));
        Assert.assertEquals(1, twoFive.compareTo(twoOne));

        GenomeLoc oneOne = GenomeLocParser.createGenomeLoc("chr1", 5);
        Assert.assertEquals(-1, oneOne.compareTo(twoOne));
        Assert.assertEquals(1, twoOne.compareTo(oneOne));
    }



}