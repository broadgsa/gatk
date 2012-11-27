package org.broadinstitute.sting.utils;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
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
        assertEquals(0, "chr1".compareTo(genomeLocParser.getContigInfo("chr1").getSequenceName())); // should be in the reference
    }

    @Test(expectedExceptions=ReviewedStingException.class)
    public void testParseBadString() {
        genomeLocParser.parseGenomeLoc("Bad:0-1");
    }

    @Test
    public void testContigHasColon() {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(net.sf.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        SAMSequenceRecord rec = new SAMSequenceRecord("c:h:r1", 10);
        rec.setSequenceLength(10);
        dict.addSequence(rec);
        header.setSequenceDictionary(dict);

        final GenomeLocParser myGenomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        GenomeLoc loc = myGenomeLocParser.parseGenomeLoc("c:h:r1:4-5");
        assertEquals(0, loc.getContigIndex());
        assertEquals(loc.getStart(), 4);
        assertEquals(loc.getStop(), 5);
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

    private static class FlankingGenomeLocTestData extends TestDataProvider {
        final GenomeLocParser parser;
        final int basePairs;
        final GenomeLoc original, flankStart, flankStop;

        private FlankingGenomeLocTestData(String name, GenomeLocParser parser, int basePairs, String original, String flankStart, String flankStop) {
            super(FlankingGenomeLocTestData.class, name);
            this.parser = parser;
            this.basePairs = basePairs;
            this.original = parse(parser, original);
            this.flankStart = flankStart == null ? null : parse(parser, flankStart);
            this.flankStop = flankStop == null ? null : parse(parser, flankStop);
        }

        private static GenomeLoc parse(GenomeLocParser parser, String str) {
            return "unmapped".equals(str) ? GenomeLoc.UNMAPPED : parser.parseGenomeLoc(str);
        }
    }

    @DataProvider(name = "flankingGenomeLocs")
    public Object[][] getFlankingGenomeLocs() {
        int contigLength = 10000;
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, contigLength);
        GenomeLocParser parser = new GenomeLocParser(header.getSequenceDictionary());

        new FlankingGenomeLocTestData("atStartBase1", parser, 1,
                "chr1:1", null, "chr1:2");

        new FlankingGenomeLocTestData("atStartBase50", parser, 50,
                "chr1:1", null, "chr1:2-51");

        new FlankingGenomeLocTestData("atStartRange50", parser, 50,
                "chr1:1-10", null, "chr1:11-60");

        new FlankingGenomeLocTestData("atEndBase1", parser, 1,
                "chr1:" + contigLength, "chr1:" + (contigLength - 1), null);

        new FlankingGenomeLocTestData("atEndBase50", parser, 50,
                "chr1:" + contigLength, String.format("chr1:%d-%d", contigLength - 50, contigLength - 1), null);

        new FlankingGenomeLocTestData("atEndRange50", parser, 50,
                String.format("chr1:%d-%d", contigLength - 10, contigLength),
                String.format("chr1:%d-%d", contigLength - 60, contigLength - 11),
                null);

        new FlankingGenomeLocTestData("nearStartBase1", parser, 1,
                "chr1:2", "chr1:1", "chr1:3");

        new FlankingGenomeLocTestData("nearStartRange50", parser, 50,
                "chr1:21-30", "chr1:1-20", "chr1:31-80");

        new FlankingGenomeLocTestData("nearEndBase1", parser, 1,
                "chr1:" + (contigLength - 1), "chr1:" + (contigLength - 2), "chr1:" + contigLength);

        new FlankingGenomeLocTestData("nearEndRange50", parser, 50,
                String.format("chr1:%d-%d", contigLength - 30, contigLength - 21),
                String.format("chr1:%d-%d", contigLength - 80, contigLength - 31),
                String.format("chr1:%d-%d", contigLength - 20, contigLength));

        new FlankingGenomeLocTestData("beyondStartBase1", parser, 1,
                "chr1:3", "chr1:2", "chr1:4");

        new FlankingGenomeLocTestData("beyondStartRange50", parser, 50,
                "chr1:101-200", "chr1:51-100", "chr1:201-250");

        new FlankingGenomeLocTestData("beyondEndBase1", parser, 1,
                "chr1:" + (contigLength - 3),
                "chr1:" + (contigLength - 4),
                "chr1:" + (contigLength - 2));

        new FlankingGenomeLocTestData("beyondEndRange50", parser, 50,
                String.format("chr1:%d-%d", contigLength - 200, contigLength - 101),
                String.format("chr1:%d-%d", contigLength - 250, contigLength - 201),
                String.format("chr1:%d-%d", contigLength - 100, contigLength - 51));

        new FlankingGenomeLocTestData("unmapped", parser, 50,
                "unmapped", null, null);

        new FlankingGenomeLocTestData("fullContig", parser, 50,
                "chr1", null, null);

        return FlankingGenomeLocTestData.getTests(FlankingGenomeLocTestData.class);
    }

    @Test(dataProvider = "flankingGenomeLocs")
    public void testCreateGenomeLocAtStart(FlankingGenomeLocTestData data) {
        GenomeLoc actual = data.parser.createGenomeLocAtStart(data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.flankStart);
        assertEquals(actual, data.flankStart, description);
    }

    @Test(dataProvider = "flankingGenomeLocs")
    public void testCreateGenomeLocAtStop(FlankingGenomeLocTestData data) {
        GenomeLoc actual = data.parser.createGenomeLocAtStop(data.original, data.basePairs);
        String description = String.format("%n      name: %s%n  original: %s%n    actual: %s%n  expected: %s%n",
                data.toString(), data.original, actual, data.flankStop);
        assertEquals(actual, data.flankStop, description);
    }
}
