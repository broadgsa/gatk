// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.


import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
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
    public void init() throws FileNotFoundException {
        // sequence
        GenomeLocParserTestUtils.clearSequenceDictionary();
        seq = new IndexedFastaSequenceFile(new File(hg18Reference));
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
        Assert.assertEquals(1, locOne.getContigIndex());
        Assert.assertEquals("chr1", locOne.getContig());

        GenomeLoc locX = GenomeLocParser.createGenomeLoc("chrX",1,1);
        Assert.assertEquals(23, locX.getContigIndex());
        Assert.assertEquals("chrX", locX.getContig());

        GenomeLoc locNumber = GenomeLocParser.createGenomeLoc(1,1,1);
        Assert.assertEquals(1, locNumber.getContigIndex());
        Assert.assertEquals("chr1", locNumber.getContig());
        Assert.assertEquals(0, locOne.compareTo(locNumber));

    }

    @Test
    public void testCompareTo() {
        logger.warn("Executing testCompareTo");
        GenomeLoc twoOne = GenomeLocParser.createGenomeLoc("chr2", 1);
        GenomeLoc twoFive = GenomeLocParser.createGenomeLoc("chr2", 5);
        GenomeLoc twoOtherFive = GenomeLocParser.createGenomeLoc("chr2", 5);
        Assert.assertEquals(twoFive.compareTo(twoOtherFive), 0);

        Assert.assertEquals(twoOne.compareTo(twoFive), -1);
        Assert.assertEquals(twoFive.compareTo(twoOne), 1);

        GenomeLoc oneOne = GenomeLocParser.createGenomeLoc("chr1", 5);
        Assert.assertEquals(oneOne.compareTo(twoOne), -1);
        Assert.assertEquals(twoOne.compareTo(oneOne), 1);
    }



}
