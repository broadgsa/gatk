package org.broadinstitute.sting.utils.fasta;

import edu.mit.broad.picard.reference.ReferenceSequence;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.BaseTest;
import org.junit.*;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 11, 2009
 * Time: 2:32:52 PM
 */
public class FastaSequenceFile2Test extends BaseTest {

    private static String sequenceFileName;
    private FastaSequenceFile2 sequenceFile = null;

    private final String firstBasesOfChrM = "GATCACAGGTCTATCACCCT";
    private final String firstBasesOfChr1 = "taaccctaaccctaacccta";
    private final String firstBasesOfChr8 = "GCAATTATGACACAAAAAAT";

    @BeforeClass
    public static void initialize() {
         sequenceFileName = seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    }

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
        logger.warn("Executing testOpenFile");

        long startTime = System.currentTimeMillis();
        Assert.assertNotNull( sequenceFile );
        long endTime = System.currentTimeMillis();

        System.err.printf("testOpenFile runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testFirstSequence() {
        logger.warn("Executing testFirstSequence");

        long startTime = System.currentTimeMillis();
        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chrM");
        Assert.assertEquals( "First n bases of chrM are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChrM.length() ),
                             firstBasesOfChrM );
        long endTime = System.currentTimeMillis();

        System.err.printf("testFirstSequence runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testNextSequence() {
        logger.warn("Executing testNextSequence");

        long startTime = System.currentTimeMillis();

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

        long endTime = System.currentTimeMillis();

        System.err.printf("testNextSequence runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testSeekToSequence() {
        logger.warn("Executing testSeekToSequence");

        long startTime = System.currentTimeMillis();

        boolean success = sequenceFile.seekToContig("chr8");
        Assert.assertTrue("Seek to seq chr8 failed", success );

        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chr8");
        Assert.assertEquals( "First n bases of chrc are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChr8.length() ),
                             firstBasesOfChr8 );

        long endTime = System.currentTimeMillis();

        System.err.printf("testSeekToSequence runtime: %dms%n", (endTime - startTime)) ;
    }

    // TODO: Is NullPointerException *really* the right exception when a sequence is missing?
    @Test(expected=NullPointerException.class)
    public void testSeekToMissingSequence() {
        logger.warn("Executing testSeekToMissingSequence");

        long startTime = 0L, endTime = 0L;

        try {
            startTime = System.currentTimeMillis();
            boolean success = sequenceFile.seekToContig("absent");
        }
        finally {
            endTime = System.currentTimeMillis();
            System.err.printf("testSeekToMissingSequence runtime: %dms%n", (endTime - startTime)) ;
        }
    }

    @Test
    public void testSeekBackward() {
        logger.warn("Executing testSeekBackward");

        long startTime = System.currentTimeMillis();

        boolean success = sequenceFile.seekToContig("chr9");
        Assert.assertTrue("Unable to seek to contig 'chr9'", success);

        success = sequenceFile.seekToContig("chr8",true);
        Assert.assertTrue("Unable to seek backward to contig 'chr8'", success);

        ReferenceSequence sequence = sequenceFile.nextSequence();
        Assert.assertEquals("First sequence contig is not correct", sequence.getName(), "chr8");
        Assert.assertEquals( "First n bases of chrc are incorrect",
                             StringUtil.bytesToString( sequence.getBases(), 0, firstBasesOfChr8.length() ),
                             firstBasesOfChr8 );

        long endTime = System.currentTimeMillis();

        System.err.printf("testSeekBackward runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testInvalidSeekBackward() {
        logger.warn("Executing testInvalidSeekBackward");

        long startTime = System.currentTimeMillis();

        boolean success = sequenceFile.seekToContig("chr9");
        Assert.assertTrue("Unable to seek to contig 'chr9'", success);

        success = sequenceFile.seekToContig("chr8");
        Assert.assertFalse("Unable to seek backward to contig 'chr8'", success);

        long endTime = System.currentTimeMillis();

        System.err.printf("testInvalidSeekBackward runtime: %dms%n", (endTime - startTime)) ;
    }

    @Test
    public void testSimultaneousAccess() {
        logger.warn("Executing testSimultaneousAccess");

        long startTime = System.currentTimeMillis();

        //        FastaSequenceFile2 other = (FastaSequenceFile2)sequenceFile.clone();

        sequenceFile.seekToContig("chr1");
        ReferenceSequence chr1 = sequenceFile.nextSequence();

//        other.seekToContig("chr8");
//        ReferenceSequence chr8 = other.nextSequence();

//        System.err.printf( "sequenceFile contig: %s%n", sequenceFile.getContigName() );
//        System.err.printf( "other contig: %s%n", other.getContigName() );

        long endTime = System.currentTimeMillis();

        System.err.printf("testSimultaneousAccess runtime: %dms%n", (endTime - startTime)) ;
    }
}
