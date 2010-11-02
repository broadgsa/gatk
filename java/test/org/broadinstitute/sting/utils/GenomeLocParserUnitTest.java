package org.broadinstitute.sting.utils;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
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
    @Test(expectedExceptions=ReviewedStingException.class)
    public void testUnsetupException() {
        SAMSequenceDictionary contigInfoCache = GenomeLocParser.contigInfo;
        GenomeLocParser.contigInfo = null;
        try {
            GenomeLocParser.createGenomeLoc(0, 0, 0);
        }
        finally {
            GenomeLocParser.contigInfo = contigInfoCache;
        }
    }

    @BeforeClass
    public void init() {
        GenomeLocParserTestUtils.clearSequenceDictionary();
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
    }

    @Test
    public void testKnownContigOrder() {
        if(true)
            throw new ReviewedStingException("Forced fail to test Bamboo pipeline");        
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        GenomeLocParser.contigInfo = null;
        // assert that it's false when the contig ordering is not setup
        assertTrue(!GenomeLocParser.hasKnownContigOrdering());
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        // assert that it's true when it is setup
        assertTrue(GenomeLocParser.hasKnownContigOrdering());
    }

    @Test(expectedExceptions=RuntimeException.class)
    public void testGetContigIndex() {
        assertEquals(GenomeLocParser.getContigIndex("blah",true), -1); // should not be in the reference
    }                

    @Test
    public void testGetContigIndexValid() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(GenomeLocParser.getContigIndex("chr1",true), 0); // should be in the reference
    }

    @Test
    public void testGetContigInfoUnknownContig() {
        assertEquals(null, GenomeLocParser.getContigInfo("blah")); // should be in the reference
    }


    @Test
    public void testGetContigInfoKnownContig() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(0, "chr1".compareTo(GenomeLocParser.getContigInfo("chr1").getSequenceName())); // should be in the reference
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testParseBadString() {
        GenomeLocParser.parseGenomeLoc("Bad:0-1");
    }

    @Test
    public void testParseGoodString() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-100");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }
    
    @Test
    public void testCreateGenomeLoc1() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1, 100);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc1point5() { // in honor of VAAL!
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc2() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1, 100);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 100);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc3() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc4() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1);
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStop(), 1);
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testCreateGenomeLoc5() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(0, 1, 100);
        GenomeLoc copy = GenomeLocParser.createGenomeLoc(loc);
        assertEquals(0, copy.getContigIndex());
        assertEquals(copy.getStop(), 100);
        assertEquals(copy.getStart(), 1);
    }

    @Test
    public void testGenomeLocPlusSign() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1+");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test
    public void testGenomeLocParseOnlyChrome() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocParseOnlyBadChrome() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr12");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocBad() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocBad2() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1-500-0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testGenomeLocBad3() {
        GenomeLoc loc = GenomeLocParser.parseGenomeLoc("chr1:1--0");
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(loc.getStop(), 10); // the size
        assertEquals(loc.getStart(), 1);
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
