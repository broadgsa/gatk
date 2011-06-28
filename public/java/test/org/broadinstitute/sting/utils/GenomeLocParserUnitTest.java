package org.broadinstitute.sting.utils;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

/**
 * @author aaron
 *         <p/>
 *         Class GenomeLocParserUnitTest
 *         <p/>
 *         Test out the functionality of the new genome loc parser
 */
public class GenomeLocParserUnitTest extends BaseTest {
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @Test(expectedExceptions=UserException.MalformedGenomeLoc.class)
    public void testGetContigIndex() {
        assertEquals(genomeLocParser.getContigIndex("blah"), -1); // should not be in the reference
    }                

    @Test
    public void testGetContigIndexValid() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(genomeLocParser.getContigIndex("chr1"), 0); // should be in the reference
    }

    @Test(expectedExceptions=UserException.class)
    public void testGetContigInfoUnknownContig1() {
        assertEquals(null, genomeLocParser.getContigInfo("blah")); // should *not* be in the reference
    }

    @Test(expectedExceptions=UserException.class)
    public void testGetContigInfoUnknownContig2() {
        assertEquals(null, genomeLocParser.getContigInfo(null)); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoUnknownContig1() {
        assertEquals(false, genomeLocParser.contigIsInDictionary("blah")); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoUnknownContig2() {
        assertEquals(false, genomeLocParser.contigIsInDictionary(null)); // should *not* be in the reference
    }

    @Test()
    public void testHasContigInfoKnownContig() {
        assertEquals(true, genomeLocParser.contigIsInDictionary("chr1")); // should be in the reference
    }

    @Test
    public void testGetContigInfoKnownContig() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(0, "chr1".compareTo(genomeLocParser.getContigInfo("chr1").getSequenceName())); // should be in the reference
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testParseBadString() {
        genomeLocParser.parseGenomeLoc("Bad:0-1");
    }

    @Test
    public void testParseGoodString() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1-10");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 10);
        assertEquals(loc.getStart(), 1);
    }
    
    @Test
    public void testCreateGenomeLoc1() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 1, 100);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc1point5() { // in honor of VAAL!
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc2() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 1, 100);
        assertEquals("chr1", loc.getContig());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc3() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 1);
        assertEquals("chr1", loc.getContig());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc4() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 1);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc5() {
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 1, 100);
        GenomeLoc copy = genomeLocParser.createGenomeLoc(loc.getContig(),loc.getStart(),loc.getStop());
        assertEquals(0, copy.getContigIndex());
        assertEquals(copy.getStop(), 100);
        assertEquals(copy.getStart(), 1);
    }

    @Test
    public void testGenomeLocPlusSign() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1+");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testGenomeLocParseOnlyChrome() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocParseOnlyBadChrome() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr12");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocBad() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1-");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.class)
    public void testGenomeLocBad2() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1-500-0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=UserException.class)
    public void testGenomeLocBad3() {
        GenomeLoc loc = genomeLocParser.parseGenomeLoc("chr1:1--0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    // test out the validating methods
    @Test
    public void testValidationOfGenomeLocs() {
        assertTrue(genomeLocParser.isValidGenomeLoc("chr1",1,1));
        assertTrue(!genomeLocParser.isValidGenomeLoc("chr2",1,1)); // shouldn't have an entry
        assertTrue(!genomeLocParser.isValidGenomeLoc("chr1",1,11)); // past the end of the contig
        assertTrue(!genomeLocParser.isValidGenomeLoc("chr1",-1,10)); // bad start
        assertTrue(!genomeLocParser.isValidGenomeLoc("chr1",1,-2)); // bad stop
        assertTrue(!genomeLocParser.isValidGenomeLoc("chr1",10,11)); // bad start, past end
    }
}
