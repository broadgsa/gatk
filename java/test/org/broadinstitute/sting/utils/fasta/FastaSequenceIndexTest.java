package org.broadinstitute.sting.utils.fasta;

import org.junit.BeforeClass;
import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;
import org.apache.log4j.BasicConfigurator;
import org.broadinstitute.sting.BaseTest;

import java.io.File;
import java.io.FileNotFoundException;

import edu.mit.broad.picard.PicardException;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 14, 2009
 * Time: 10:34:15 AM
 * To change this template use File | Settings | File Templates.
 */
public class FastaSequenceIndexTest extends BaseTest {
    private static String sequenceIndexName = null;
    private FastaSequenceIndex sequenceIndex = null;

    @BeforeClass
    public static void initialize() {
        sequenceIndexName = seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai";
    }

    @Before
    public void doForEachTest() throws FileNotFoundException {
        sequenceIndex = new FastaSequenceIndex( new File(sequenceIndexName) );
    }

    @Test
    public void testInitialContig() {
        Assert.assertTrue("Contig chrM is not present", sequenceIndex.hasIndexEntry("chrM"));
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chrM");
        Assert.assertEquals("Contig chrM name is incorrect",entry.getContig(),"chrM");
        Assert.assertEquals("Contig chrM location is incorrect",entry.getLocation(),6L);
        Assert.assertEquals("Contig chrM size is incorrect",entry.getSize(),16571L);
        Assert.assertEquals("Contig chrM bases per line is incorrect",entry.getBasesPerLine(),50);
        Assert.assertEquals("Contig chrM bytes per line is incorrect",entry.getBytesPerLine(),51);
    }

    @Test
    public void testMiddleContig() {
        Assert.assertTrue("Contig chr8 is not present", sequenceIndex.hasIndexEntry("chr8"));
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chr8");
        Assert.assertEquals("Contig chr8 name is incorrect",entry.getContig(),"chr8");
        Assert.assertEquals("Contig chr8 location is incorrect",entry.getLocation(),1419403101L);
        Assert.assertEquals("Contig chr8 size is incorrect",entry.getSize(),146274826L);
        Assert.assertEquals("Contig chr8 bases per line is incorrect",entry.getBasesPerLine(),50);
        Assert.assertEquals("Contig chr8 bytes per line is incorrect",entry.getBytesPerLine(),51);
    }

    @Test
    public void testLastContig() {
        Assert.assertTrue("Contig chrX_random is not present", sequenceIndex.hasIndexEntry("chrX_random"));
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chrX_random");
        Assert.assertEquals("Contig chrX_random name is incorrect",entry.getContig(),"chrX_random");
        Assert.assertEquals("Contig chrX_random location is incorrect",entry.getLocation(),3156698441L);
        Assert.assertEquals("Contig chrX_random size is incorrect",entry.getSize(),1719168L);
        Assert.assertEquals("Contig chrX_random bases per line is incorrect",entry.getBasesPerLine(),50);
        Assert.assertEquals("Contig chrX_random bytes per line is incorrect",entry.getBytesPerLine(),51);
    }

    @Test
    public void testAllContigsPresent() {
        Assert.assertTrue("Contig chrM is not present", sequenceIndex.hasIndexEntry("chrM"));
        Assert.assertTrue("Contig chr1 is not present", sequenceIndex.hasIndexEntry("chr1"));
        Assert.assertTrue("Contig chr2 is not present", sequenceIndex.hasIndexEntry("chr2"));
        Assert.assertTrue("Contig chr3 is not present", sequenceIndex.hasIndexEntry("chr3"));
        Assert.assertTrue("Contig chr4 is not present", sequenceIndex.hasIndexEntry("chr4"));
        Assert.assertTrue("Contig chr5 is not present", sequenceIndex.hasIndexEntry("chr5"));
        Assert.assertTrue("Contig chr6 is not present", sequenceIndex.hasIndexEntry("chr6"));
        Assert.assertTrue("Contig chr7 is not present", sequenceIndex.hasIndexEntry("chr7"));
        Assert.assertTrue("Contig chr8 is not present", sequenceIndex.hasIndexEntry("chr8"));
        Assert.assertTrue("Contig chr9 is not present", sequenceIndex.hasIndexEntry("chr9"));
        Assert.assertTrue("Contig chr10 is not present", sequenceIndex.hasIndexEntry("chr10"));
        Assert.assertTrue("Contig chr11 is not present", sequenceIndex.hasIndexEntry("chr11"));
        Assert.assertTrue("Contig chr12 is not present", sequenceIndex.hasIndexEntry("chr12"));
        Assert.assertTrue("Contig chr13 is not present", sequenceIndex.hasIndexEntry("chr13"));
        Assert.assertTrue("Contig chr14 is not present", sequenceIndex.hasIndexEntry("chr14"));
        Assert.assertTrue("Contig chr15 is not present", sequenceIndex.hasIndexEntry("chr15"));
        Assert.assertTrue("Contig chr16 is not present", sequenceIndex.hasIndexEntry("chr16"));
        Assert.assertTrue("Contig chr17 is not present", sequenceIndex.hasIndexEntry("chr17"));
        Assert.assertTrue("Contig chr18 is not present", sequenceIndex.hasIndexEntry("chr18"));
        Assert.assertTrue("Contig chr19 is not present", sequenceIndex.hasIndexEntry("chr19"));
        Assert.assertTrue("Contig chr20 is not present", sequenceIndex.hasIndexEntry("chr20"));
        Assert.assertTrue("Contig chr21 is not present", sequenceIndex.hasIndexEntry("chr21"));
        Assert.assertTrue("Contig chr22 is not present", sequenceIndex.hasIndexEntry("chr22"));
        Assert.assertTrue("Contig chrX is not present", sequenceIndex.hasIndexEntry("chrX"));
        Assert.assertTrue("Contig chrY is not present", sequenceIndex.hasIndexEntry("chrY"));
        Assert.assertTrue("Contig chr1_random is not present", sequenceIndex.hasIndexEntry("chr1_random"));
        Assert.assertTrue("Contig chr2_random is not present", sequenceIndex.hasIndexEntry("chr2_random"));
        Assert.assertTrue("Contig chr3_random is not present", sequenceIndex.hasIndexEntry("chr3_random"));
        Assert.assertTrue("Contig chr4_random is not present", sequenceIndex.hasIndexEntry("chr4_random"));
        Assert.assertTrue("Contig chr5_random is not present", sequenceIndex.hasIndexEntry("chr5_random"));
        Assert.assertTrue("Contig chr6_random is not present", sequenceIndex.hasIndexEntry("chr6_random"));
        Assert.assertTrue("Contig chr7_random is not present", sequenceIndex.hasIndexEntry("chr7_random"));
        Assert.assertTrue("Contig chr8_random is not present", sequenceIndex.hasIndexEntry("chr8_random"));
        Assert.assertTrue("Contig chr9_random is not present", sequenceIndex.hasIndexEntry("chr9_random"));
        Assert.assertTrue("Contig chr10_random is not present", sequenceIndex.hasIndexEntry("chr10_random"));
        Assert.assertTrue("Contig chr11_random is not present", sequenceIndex.hasIndexEntry("chr11_random"));
        Assert.assertTrue("Contig chr13_random is not present", sequenceIndex.hasIndexEntry("chr13_random"));
        Assert.assertTrue("Contig chr15_random is not present", sequenceIndex.hasIndexEntry("chr15_random"));
        Assert.assertTrue("Contig chr16_random is not present", sequenceIndex.hasIndexEntry("chr16_random"));
        Assert.assertTrue("Contig chr17_random is not present", sequenceIndex.hasIndexEntry("chr17_random"));
        Assert.assertTrue("Contig chr18_random is not present", sequenceIndex.hasIndexEntry("chr18_random"));
        Assert.assertTrue("Contig chr19_random is not present", sequenceIndex.hasIndexEntry("chr19_random"));
        Assert.assertTrue("Contig chr21_random is not present", sequenceIndex.hasIndexEntry("chr21_random"));
        Assert.assertTrue("Contig chr22_random is not present", sequenceIndex.hasIndexEntry("chr22_random"));
        Assert.assertTrue("Contig chrX_random is not present", sequenceIndex.hasIndexEntry("chrX_random"));
    }

    @Test
    public void testHasInvalidEntry() {
        Assert.assertFalse("Found an invalid entry", sequenceIndex.hasIndexEntry("invalid"));
    }

    @Test(expected= PicardException.class)
    public void testGetInvalidEntry() {
        sequenceIndex.getIndexEntry("invalid");
    }

}
