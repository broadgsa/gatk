package org.broadinstitute.sting.utils;

import org.junit.Test;
import org.junit.Assert;
import org.junit.Before;
import org.junit.After;

import java.io.File;

import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 11, 2009
 * Time: 2:32:52 PM
 */
public class FastaSequenceFile2Test {

    private final String sequenceFileName = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    private FastaSequenceFile2 sequenceFile = null;

    private final String firstBasesOfChrM = "GATCACAGGTCTATCACCCT";
    private final String firstBasesOfChr1 = "taaccctaaccctaacccta";
    private final String firstBasesOfChr8 = "GCAATTATGACACAAAAAAT";

    @Before
    public void doForEachTest() {
        sequenceFile = new FastaSequenceFile2( new File(sequenceFileName) );
    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @After
    public void undoForEachTest() {
        sequenceFile = null;
    }


    @Test
    public void testOpenFile() {
        Assert.assertNotNull( sequenceFile );
    }

    @Test
    public void testFirstSequence() {
        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals( "First n bases of chrM are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChrM.length() ),
                             firstBasesOfChrM );
    }

    @Test
    public void testNextSequence() {
        ReferenceSequence sequence = null;

        // Advance to chrM.
        sequence = sequenceFile.nextSequence();
        sequence = sequenceFile.nextSequence();

        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chr1");

        // Workaround: bytesToString for chr1 of the fasta file we've picked doesn't appear to work.
        // TODO: Report this as sam-jdk bug.
        byte[] firstOfChr1 = StringUtil.stringToBytes(firstBasesOfChr1);
        byte[] firstOfSequence = new byte[firstBasesOfChr1.length()];
        System.arraycopy(sequence.getBases(), 0, firstOfSequence, 0, firstOfSequence.length );

        Assert.assertArrayEquals("First bases of chr1 are not correct", firstOfChr1, firstOfSequence );
    }

    @Test
    public void testSeekToSequence() {
        boolean success = sequenceFile.seekToContig("chr8");
        Assert.assertTrue("Seek to seq chr8 failed", success );

        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chr8");
        Assert.assertEquals( "First n bases of chrc are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChr8.length() ),
                             firstBasesOfChr8 );
    }

    // TODO: Is NullPointerException *really* the right exception when a sequence is missing?
    @Test(expected=NullPointerException.class)
    public void testSeekToMissingSequence() {
        boolean success = sequenceFile.seekToContig("absent");
    }

    @Test
    public void testSeekBackward() {
        boolean success = sequenceFile.seekToContig("chr9");
        Assert.assertTrue("Unable to seek to contig 'chr9'", success);

        success = sequenceFile.seekToContig("chr8",true);
        Assert.assertTrue("Unable to seek backward to contig 'chr8'", success);

        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chr8");
        Assert.assertEquals( "First n bases of chrc are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChr8.length() ),
                             firstBasesOfChr8 );        
    }

    @Test
    public void testInvalidSeekBackward() {
        boolean success = sequenceFile.seekToContig("chr9");
        Assert.assertTrue("Unable to seek to contig 'chr9'", success);

        success = sequenceFile.seekToContig("chr8");
        Assert.assertFalse("Unable to seek backward to contig 'chr8'", success);
    }
}
