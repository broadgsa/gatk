package org.broadinstitute.sting.utils;

import static junit.framework.Assert.assertTrue;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class GenomeLocParserTest
 *         <p/>
 *         Test out the functionality of the new genome loc parser
 */
public class GenomeLocParserTest extends BaseTest {

    @Test(expected = StingException.class)
    public void testUnsetupException() {
        GenomeLocParser.createGenomeLoc(0, 0, 0);
    }

    @Test
    public void testKnownContigOrder() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        // assert that it's false when the contig ordering is not setup
        assertTrue(!GenomeLocParser.hasKnownContigOrdering());
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        // assert that it's true when it is setup
        assertTrue(GenomeLocParser.hasKnownContigOrdering());
    }

    @Test(expected = RuntimeException.class)
    public void testGetContigIndex() {
        assertEquals(-1, GenomeLocParser.getContigIndex("blah")); // should be in the reference
    }

    @Test
    public void testGetContigIndexValid() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 10);
        assertEquals(0, GenomeLocParser.getContigIndex("chr1")); // should be in the reference
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

    @Test(expected = RuntimeException.class)
    public void testParseBadLocations() {
        GenomeLocParser.parseGenomeLocs("chr1:1-1;badChr:1-0");
    }

    @Test
    public void testParseGoodLocations() {
        GenomeLocParser.parseGenomeLocs("chr1:1-1;chr1:5-9");
    }

    @Test(expected = RuntimeException.class)
    public void testParseGoodLocationsTooManySemiColons() {
        GenomeLocParser.parseGenomeLocs("chr1:1-1;;chr1:5-9;");
    }

    @Test
    public void testCreateGenomeLoc1() {
        GenomeLoc loc = GenomeLocParser.createGenomeLoc("chr1", 1, 100);
        assertEquals(loc.getContigIndex(), 0);
        assertEquals(100, loc.getStop());
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
    public void testGenomeLocParserList() {
        long start = System.currentTimeMillis();
        List<GenomeLoc> parsedIntervals = GenomeAnalysisEngine.parseIntervalRegion(Arrays.asList(new String[]{"/humgen/gsa-scr1/GATK_Data/Validation_Data/bigChr1IntervalList.list"}));
        Collections.sort(parsedIntervals);
        LinkedList<GenomeLoc> loc = new LinkedList<GenomeLoc>(GenomeLocParser.mergeOverlappingLocations(parsedIntervals));
        long stop = System.currentTimeMillis();
        logger.warn("Elapsed time = " + (stop - start));
    }
}
