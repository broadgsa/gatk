// our package
package org.broadinstitute.sting.utils;


// the imports for unit testing.


import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;

/**
 * Basic unit test for GenomeLoc
 */
public class GenomeLocUnitTest extends BaseTest {
    private static ReferenceSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testIsBetween() {
        logger.warn("Executing testIsBetween");

        GenomeLoc locMiddle = genomeLocParser.createGenomeLoc("chr1", 3, 3);

        GenomeLoc locLeft = genomeLocParser.createGenomeLoc("chr1", 1, 1);
        GenomeLoc locRight = genomeLocParser.createGenomeLoc("chr1", 5, 5);

        Assert.assertTrue(locMiddle.isBetween(locLeft, locRight));
        Assert.assertFalse(locLeft.isBetween(locMiddle, locRight));
        Assert.assertFalse(locRight.isBetween(locLeft, locMiddle));

    }
    @Test
    public void testContigIndex() {
        logger.warn("Executing testContigIndex");
        GenomeLoc locOne = genomeLocParser.createGenomeLoc("chr1",1,1);
        Assert.assertEquals(1, locOne.getContigIndex());
        Assert.assertEquals("chr1", locOne.getContig());

        GenomeLoc locX = genomeLocParser.createGenomeLoc("chrX",1,1);
        Assert.assertEquals(23, locX.getContigIndex());
        Assert.assertEquals("chrX", locX.getContig());

        GenomeLoc locNumber = genomeLocParser.createGenomeLoc(seq.getSequenceDictionary().getSequence(1).getSequenceName(),1,1);
        Assert.assertEquals(1, locNumber.getContigIndex());
        Assert.assertEquals("chr1", locNumber.getContig());
        Assert.assertEquals(0, locOne.compareTo(locNumber));

    }

    @Test
    public void testCompareTo() {
        logger.warn("Executing testCompareTo");
        GenomeLoc twoOne = genomeLocParser.createGenomeLoc("chr2", 1);
        GenomeLoc twoFive = genomeLocParser.createGenomeLoc("chr2", 5);
        GenomeLoc twoOtherFive = genomeLocParser.createGenomeLoc("chr2", 5);
        Assert.assertEquals(twoFive.compareTo(twoOtherFive), 0);

        Assert.assertEquals(twoOne.compareTo(twoFive), -1);
        Assert.assertEquals(twoFive.compareTo(twoOne), 1);

        GenomeLoc oneOne = genomeLocParser.createGenomeLoc("chr1", 5);
        Assert.assertEquals(oneOne.compareTo(twoOne), -1);
        Assert.assertEquals(twoOne.compareTo(oneOne), 1);
    }



}
