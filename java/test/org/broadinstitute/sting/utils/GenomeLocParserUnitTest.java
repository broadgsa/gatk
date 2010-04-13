package org.broadinstitute.sting.utils;


import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.BeforeClass;
import org.junit.Test;

/**
 * @author aaron
 *         <p/>
 *         Class GenomeLocParserUnitTest
 *         <p/>
 *         Test out the functionality of the new genome loc parser
 */
public class GenomeLocParserUnitTest extends BaseTest {
    @Test(expected = StingException.class)
    public void testUnsetupException() {
        GenomeLocParser.contigInfo = null;
        GenomeLocParser.createGenomeLoc(0, 0, 0);
    }

    @BeforeClass
    public static void init() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Test
    public void testKnownContigOrder() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        GenomeLocParser.contigInfo = null;
        // assert that it's false when the contig ordering is not setup
        assertTrue(!GenomeLocParser.hasKnownContigOrdering());
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        // assert that it's true when it is setup
        assertTrue(GenomeLocParser.hasKnownContigOrdering());
    }

    @Test(expected = RuntimeException.class)
    public void testGetContigIndex() {
        assertEquals(-1, GenomeLocParser.getContigIndex("blah",true)); // should not be in the reference
    }                

    @Test
    public void testGetContigIndexValid() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(0, GenomeLocParser.getContigIndex("chr1",true)); // should be in the reference
    }

    @Test
    public void testGetContigInfoUnknownContig() {
        assertEquals(null, GenomeLocParser.getContigInfo("blah")); // should be in the reference
    }


    @Test
    public void testGetContigInfoKnownContig() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals("chr1".compareTo(GenomeLocParser.getContigInfo("chr1").getSequenceName()), 0); // should be in the reference
    }

    @Test(expected = StingException.class)
    public void testParseBadString() {
        GenomeLocParser.parseGenomeLoc("Bad:0-1");
    }

    @Test
    public void testParseGoodString() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-100");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(100, loc.getStop());
        assertEquals(1, loc.getStart());
    }
    
    @Test
    public void testCreateGenomeLoc1() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1, 100);
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(100, loc.getStop());
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testCreateGenomeLoc1point5() { // in honor of VAAL!
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(1, loc.getStop());
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testCreateGenomeLoc2() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1, 100);
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(100, loc.getStop());
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testCreateGenomeLoc3() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1);
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(1, loc.getStop());
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testCreateGenomeLoc4() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1);
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(1, loc.getStop());
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testCreateGenomeLoc5() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1, 100);
        GenomeLoc copy = GenomeLocParser.createGenomeLoc(loc);
        assertEquals(copy.getContigIndex(), 0);
        assertEquals(100, copy.getStop());
        assertEquals(1, copy.getStart());
    }

    @Test
    public void testGenomeLocPlusSign() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1+");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    @Test
    public void testGenomeLocParseOnlyChrome() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    @Test(expected = StingException.class)
    public void testGenomeLocParseOnlyBadChrome() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr12");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    @Test(expected = StingException.class)
    public void testGenomeLocBad() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    @Test(expected = StingException.class)
    public void testGenomeLocBad2() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-500-0");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    @Test(expected = StingException.class)
    public void testGenomeLocBad3() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1--0");
        assertEquals(0, loc.getContigIndex());
        assertEquals(10, loc.getStop()); // the size
        assertEquals(1, loc.getStart());
    }

    // test out the validating methods
    @Test
    public void testValidationOfGenomeLocs() {
        assertTrue(GenomeLocParser.validGenomeLoc("chr1",1,1));
        assertTrue(!GenomeLocParser.validGenomeLoc("chr2",1,1)); // shouldn't have an entry
        assertTrue(!GenomeLocParser.validGenomeLoc("chr1",1,11)); // past the end of the contig
        assertTrue(!GenomeLocParser.validGenomeLoc("chr1",-1,10)); // bad start
        assertTrue(!GenomeLocParser.validGenomeLoc("chr1",1,-2)); // bad stop
        assertTrue(!GenomeLocParser.validGenomeLoc("chr1",10,11)); // bad start, past end

        assertTrue(GenomeLocParser.validGenomeLoc(0,1,1));
        assertTrue(!GenomeLocParser.validGenomeLoc(1,1,1)); // shouldn't have an entry
        assertTrue(!GenomeLocParser.validGenomeLoc(0,1,11)); // past the end of the contig
        assertTrue(!GenomeLocParser.validGenomeLoc(-1,0,10)); // bad start
        assertTrue(!GenomeLocParser.validGenomeLoc(0,1,-2)); // bad stop
        assertTrue(!GenomeLocParser.validGenomeLoc(0,10,11)); // bad start, past end

    }
}
