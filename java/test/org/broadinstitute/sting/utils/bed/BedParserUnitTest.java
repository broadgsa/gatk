package org.broadinstitute.sting.utils.bed;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


import java.io.File;
import java.util.List;

import net.sf.picard.reference.IndexedFastaSequenceFile;


public class BedParserUnitTest extends BaseTest {

    private static IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;
    private File bedFile = new File("testdata/sampleBedFile.bed");

    @BeforeClass
    public void beforeTests() {
        seq = new IndexedFastaSequenceFile(new File(b36KGReference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test
    public void testLoadBedFile() {
        BedParser parser = new BedParser(genomeLocParser,bedFile);
        List<GenomeLoc> location = parser.getLocations();
        Assert.assertEquals(location.size(), 4);
    }

    @Test
    public void testBedParsing() {
        BedParser parser = new BedParser(genomeLocParser,bedFile);
        List<GenomeLoc> location = parser.getLocations();
        Assert.assertEquals(location.size(), 4);
        Assert.assertTrue(location.get(0).getContig().equals("20"));
        Assert.assertTrue(location.get(1).getContig().equals("20"));
        Assert.assertTrue(location.get(2).getContig().equals("22"));
        Assert.assertTrue(location.get(3).getContig().equals("22"));

        // now check the the start positions
        Assert.assertEquals(location.get(0).getStart(), 1);
        Assert.assertEquals(location.get(1).getStart(), 1002);
        Assert.assertEquals(location.get(2).getStart(), 1001);
        Assert.assertEquals(location.get(3).getStart(), 2001);

        // now check the the stop positions
        Assert.assertEquals(location.get(0).getStop(), 999);
        Assert.assertEquals(location.get(1).getStop(), 2000);
        Assert.assertEquals(location.get(2).getStop(), 5000);
        Assert.assertEquals(location.get(3).getStop(), 6000);
    }
}
