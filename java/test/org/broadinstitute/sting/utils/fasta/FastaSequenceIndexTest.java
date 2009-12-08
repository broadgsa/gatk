package org.broadinstitute.sting.utils.fasta;

import net.sf.picard.PicardException;
import org.broadinstitute.sting.BaseTest;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

/**
 * Test the fasta sequence index reader.
 */
public class FastaSequenceIndexTest extends BaseTest {
    // our basic human 18 fai
    private static String sequenceIndexName = null;
    private FastaSequenceIndex sequenceIndex = null;

    // a custom index that tests the colon, and semi-colon, and other random characters
    private static String sequenceIndexColonSemiColonTestName = null;
    private FastaSequenceIndex sequenceIndexColonSemiColonTest = null;


    @BeforeClass
    public static void initialize() {
        sequenceIndexName = seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.fai";
        sequenceIndexColonSemiColonTestName = "/humgen/gsa-scr1/GATK_Data/Validation_Data/testing.fai";
    }

    @Before
    public void doForEachTest() throws FileNotFoundException {
        sequenceIndex = new FastaSequenceIndex( new File(sequenceIndexName) );
        sequenceIndexColonSemiColonTest = new FastaSequenceIndex( new File(sequenceIndexColonSemiColonTestName) );
    }

    @Test
    public void testInitialContig() {
        logger.warn("Executing testInitialContig");

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
        logger.warn("Executing testMiddleContig");

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
        logger.warn("Executing testLastContig");

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
        logger.warn("Executing testAllContigsPresent");

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
        logger.warn("Executing testHasInvalidEntry");

        Assert.assertFalse("Found an invalid entry", sequenceIndex.hasIndexEntry("invalid"));
    }

    @Test(expected= PicardException.class)
    public void testGetInvalidEntry() {
        logger.warn("Executing testGetInvalidEntry");

        sequenceIndex.getIndexEntry("invalid");
    }

    @Test
    public void testIteration() {
        logger.warn("Executing testIteration");        

        Iterator<FastaSequenceIndexEntry> sequenceIndexEntries = sequenceIndex.iterator();

        Assert.assertEquals("Contig chrM is not present", "chrM", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr1 is not present", "chr1", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr2 is not present", "chr2", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr3 is not present", "chr3", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr4 is not present", "chr4", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr5 is not present", "chr5", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr6 is not present", "chr6", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr7 is not present", "chr7", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr8 is not present", "chr8", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr9 is not present", "chr9", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr10 is not present", "chr10", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr11 is not present", "chr11", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr12 is not present", "chr12", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr13 is not present", "chr13", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr14 is not present", "chr14", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr15 is not present", "chr15", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr16 is not present", "chr16", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr17 is not present", "chr17", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr18 is not present", "chr18", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr19 is not present", "chr19", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr20 is not present", "chr20", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr21 is not present", "chr21", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr22 is not present", "chr22", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chrX is not present", "chrX", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chrY is not present", "chrY", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr1_random is not present", "chr1_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr2_random is not present", "chr2_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr3_random is not present", "chr3_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr4_random is not present", "chr4_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr5_random is not present", "chr5_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr6_random is not present", "chr6_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr7_random is not present", "chr7_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr8_random is not present", "chr8_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr9_random is not present", "chr9_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr10_random is not present", "chr10_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr11_random is not present", "chr11_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr13_random is not present", "chr13_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr15_random is not present", "chr15_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr16_random is not present", "chr16_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr17_random is not present", "chr17_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr18_random is not present", "chr18_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr19_random is not present", "chr19_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr21_random is not present", "chr21_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chr22_random is not present", "chr22_random", sequenceIndexEntries.next().getContig());
        Assert.assertEquals("Contig chrX_random is not present", "chrX_random", sequenceIndexEntries.next().getContig());
        Assert.assertFalse("Iterator still has more entries", sequenceIndexEntries.hasNext());
    }

    @Test
    public void testSpecialCharacters() {
        /* file contents:
        chrM	16571	6	50	51
        chr1;boat	247249719	16915	50	51
        chr2:money	242951149	252211635	50	51
        chr3::;	199501827	500021813	50	51
        ;;;;;;  1234            1234            1234    1234
        file:gi|17981852|ref|NC_001807.4|    16571   2911876801      70      71
        */
        Iterator<FastaSequenceIndexEntry> sequenceIndexEntries = sequenceIndexColonSemiColonTest.iterator();
        FastaSequenceIndexEntry ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig chrM is not present","chrM",ent.getContig());
        Assert.assertEquals("Contig chrM size is not correct",16571,ent.getSize());
        Assert.assertEquals("Contig chrM location is not correct",6,ent.getLocation());
        Assert.assertEquals("Contig chrM bases per line is not correct",50,ent.getBasesPerLine());
        Assert.assertEquals("Contig chrM bytes per line is not correct",51,ent.getBytesPerLine());

        ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig chr1;boat is not present","chr1;boat",ent.getContig());
        Assert.assertEquals("Contig chr1;boat size is not correct",247249719,ent.getSize());
        Assert.assertEquals("Contig chr1;boat location is not correct",16915,ent.getLocation());
        Assert.assertEquals("Contig chr1;boat bases per line is not correct",50,ent.getBasesPerLine());
        Assert.assertEquals("Contig chr1;boat bytes per line is not correct",51,ent.getBytesPerLine());

        ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig chr2:money is not present","chr2:money",ent.getContig());
        Assert.assertEquals("Contig chr2:money size is not correct",242951149,ent.getSize());
        Assert.assertEquals("Contig chr2:money location is not correct",252211635,ent.getLocation());
        Assert.assertEquals("Contig chr2:money bases per line is not correct",50,ent.getBasesPerLine());
        Assert.assertEquals("Contig chr2:money bytes per line is not correct",51,ent.getBytesPerLine());

        ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig chr3::; is not present","chr3::;",ent.getContig());
        Assert.assertEquals("Contig chr3::; size is not correct",199501827,ent.getSize());
        Assert.assertEquals("Contig chrM location is not correct",500021813,ent.getLocation());
        Assert.assertEquals("Contig chr3::; bases per line is not correct",50,ent.getBasesPerLine());
        Assert.assertEquals("Contig chr3::; bytes per line is not correct",51,ent.getBytesPerLine());

        ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig ;;;;;;;; is not present",";;;;;;;;",ent.getContig());
        Assert.assertEquals("Contig ;;;;;;;; size is not correct",123,ent.getSize());
        Assert.assertEquals("Contig ;;;;;;;; location is not correct",234,ent.getLocation());
        Assert.assertEquals("Contig ;;;;;;;; bases per line is not correct",456,ent.getBasesPerLine());
        Assert.assertEquals("Contig ;;;;;;;; bytes per line is not correct",789,ent.getBytesPerLine());

        ent = sequenceIndexEntries.next();
        Assert.assertEquals("Contig file:gi|17981852|ref|NC_001807.4| is not present","file:gi|17981852|ref|NC_001807.4|",ent.getContig());
        Assert.assertEquals("Contig file:gi|17981852|ref|NC_001807.4| size is not correct",16571,ent.getSize());
        Assert.assertEquals("Contig file:gi|17981852|ref|NC_001807.4| location is not correct",2911876801L,ent.getLocation());
        Assert.assertEquals("Contig file:gi|17981852|ref|NC_001807.4| bases per line is not correct",70,ent.getBasesPerLine());
        Assert.assertEquals("Contig file:gi|17981852|ref|NC_001807.4| bytes per line is not correct",71,ent.getBytesPerLine());
    }
}
