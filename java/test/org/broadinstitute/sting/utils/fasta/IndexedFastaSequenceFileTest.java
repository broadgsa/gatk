package org.broadinstitute.sting.utils.fasta;

import org.junit.BeforeClass;
import org.junit.Before;
import org.junit.Test;
import org.junit.Assert;
import org.broadinstitute.sting.BaseTest;

import java.io.File;
import java.io.FileNotFoundException;

import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 14, 2009
 * Time: 2:37:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class IndexedFastaSequenceFileTest extends BaseTest {
    private static String sequenceFileName;
    private IndexedFastaSequenceFile sequenceFile = null;

    private final String firstBasesOfChrM = "GATCACAGGTCTATCACCCT";
    private final String extendedBasesOfChrM = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT" +
                                               "TTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG" +
                                               "GAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATT";    
    private final String firstBasesOfChr1 = "taaccctaaccctaacccta";
    private final String firstBasesOfChr8 = "GCAATTATGACACAAAAAAT";

    @BeforeClass
    public static void initialize() {
         sequenceFileName = seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    }

    @Before
    public void doForEachTest() throws FileNotFoundException {
        sequenceFile = new IndexedFastaSequenceFile( new File(sequenceFileName) );
    }

    @Test
    public void testOpenFile() {
        long startTime = System.currentTimeMillis();
        Assert.assertNotNull( sequenceFile );
        long endTime = System.currentTimeMillis();

        System.err.printf("testOpenFile runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testFirstSequence() {
        long startTime = System.currentTimeMillis();
        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals( "First n bases of chrM are incorrect",
                             firstBasesOfChrM,
                             StringUtil.bytesToString( sequence.getBases() ) );
        long endTime = System.currentTimeMillis();

        System.err.printf("testFirstSequence runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testFirstSequenceExtended() {
        long startTime = System.currentTimeMillis();
        ReferenceSequence sequence = sequenceFile.getSubsequenceAt("chrM",0,extendedBasesOfChrM.length());
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals( "First n bases of chrM are incorrect",
                             extendedBasesOfChrM.substring(0,110),
                             StringUtil.bytesToString( sequence.getBases(),0,110 ) );
        long endTime = System.currentTimeMillis();

        System.err.printf("testFirstSequenceExtended runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testReadStartingInCenterOfLine() {
        final int bytesToChopOff = 5;
        String truncated = extendedBasesOfChrM.substring(bytesToChopOff);

        long startTime = System.currentTimeMillis();
        ReferenceSequence sequence = sequenceFile.getSubsequenceAt("chrM", bytesToChopOff ,truncated.length() );
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals( "First n bases of chrM are incorrect",
                             truncated,
                             StringUtil.bytesToString( sequence.getBases() ) );
        long endTime = System.currentTimeMillis();

        System.err.printf("testReadStartingInCenterOfLine runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testCompleteContigRead() {
        FastaSequenceFile2 originalSequenceFile = new FastaSequenceFile2(new File(sequenceFileName));
        ReferenceSequence expectedSequence = originalSequenceFile.nextSequence();

        long startTime = System.currentTimeMillis();
        ReferenceSequence sequence = sequenceFile.getSequence("chrM");        
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals("chrM is incorrect",
                            StringUtil.bytesToString(expectedSequence.getBases(),0,4096),
                            StringUtil.bytesToString(sequence.getBases(),0,4096) );
        long endTime = System.currentTimeMillis();

        System.err.printf("testCompleteContigRead runtime: %dms%n", (endTime - startTime)) ;        
    }


}
