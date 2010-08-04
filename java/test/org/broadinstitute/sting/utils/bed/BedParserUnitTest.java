package org.broadinstitute.sting.utils.bed;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import net.sf.picard.reference.IndexedFastaSequenceFile;


public class BedParserUnitTest extends BaseTest {

    private static IndexedFastaSequenceFile seq;
    private File bedFile = new File("testdata/sampleBedFile.bed");

    @BeforeClass
    public static void beforeTests() {
        seq = new IndexedFastaSequenceFile(new File(b36KGReference));
        GenomeLocParser.setupRefContigOrdering(seq);
    }

    @Test
    public void testLoadBedFile() {
        BedParser parser = new BedParser(bedFile);
        List<GenomeLoc> location = parser.getLocations();
        Assert.assertEquals(4, location.size());
    }

    @Test
    public void testBedParsing() {
        BedParser parser = new BedParser(bedFile);
        List<GenomeLoc> location = parser.getLocations();
        Assert.assertEquals(4, location.size());
        Assert.assertTrue(location.get(0).getContig().equals("20"));
        Assert.assertTrue(location.get(1).getContig().equals("20"));
        Assert.assertTrue(location.get(2).getContig().equals("22"));
        Assert.assertTrue(location.get(3).getContig().equals("22"));

        // now check the the start positions
        Assert.assertEquals(1, location.get(0).getStart());
        Assert.assertEquals(1002, location.get(1).getStart());
        Assert.assertEquals(1001, location.get(2).getStart());
        Assert.assertEquals(2001, location.get(3).getStart());

        // now check the the stop positions
        Assert.assertEquals(999, location.get(0).getStop());
        Assert.assertEquals(2000, location.get(1).getStop());
        Assert.assertEquals(5000, location.get(2).getStop());
        Assert.assertEquals(6000, location.get(3).getStop());
    }
}
