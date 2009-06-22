package org.broadinstitute.sting.playground.gatk.walkers.indels;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.OutputTracker;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileReader;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileWriter;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.File;
import java.util.Arrays;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import junit.framework.Assert;
/**
 * User: hanna
 * Date: Jun 11, 2009
 * Time: 9:54:05 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Test the walker that injects clean reads into the bam file.
 */

public class CleanedReadInjectorTest extends BaseTest {
    /**
     * The fasta, for comparison.
     */
    protected static IndexedFastaSequenceFile sequenceFile = null;

    /**
     * Initialize the fasta.
     */
    @BeforeClass
    public static void initialize() throws FileNotFoundException {
        sequenceFile = new IndexedFastaSequenceFile( new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta") );
        GenomeLocParser.setupRefContigOrdering(sequenceFile);
    }

    @Test
    public void testNoReads() {
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader();
        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:1-10", cleanedReads, output );

        walker.initialize();
        walker.onTraversalDone(0);

        Assert.assertEquals("Too many records in output",0,output.getRecords().size());
    }

    @Test
    public void testNoCleanedReads() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord sourceRead = ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,5);

        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader();
        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:1-10", cleanedReads, output );

        int result = runWalkerOverReads(walker,sourceRead);

        Assert.assertEquals("Result of traversal is incorrect",0,result);
        Assert.assertEquals("Incorrect number of records in output",1,output.getRecords().size());
        Assert.assertEquals("Output record is incorrect",sourceRead,output.getRecords().get(0));
    }

    @Test
    public void testOneCleanedRead() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord sourceRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);

        SAMRecord cleanedRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);
        cleanedRead.setBaseQualities(getMockBaseQualityString((byte)1,cleanedRead.getReadLength()));
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader(cleanedRead);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:1-10", cleanedReads, output );

        int result = runWalkerOverReads(walker,sourceRead);

        Assert.assertEquals("Result of traversal is incorrect",1,result);
        Assert.assertEquals("Incorrect number of records in output",1,output.getRecords().size());
        Assert.assertEquals("Output record is incorrect",cleanedRead,output.getRecords().get(0));
    }

    @Test
    public void testPartialIntervalOverlap() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord sourceRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);

        SAMRecord cleanedRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);
        cleanedRead.setBaseQualities(getMockBaseQualityString((byte)1,cleanedRead.getReadLength()));
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader(cleanedRead);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:4-12", cleanedReads, output );

        int result = runWalkerOverReads(walker,sourceRead);

        Assert.assertEquals("Result of traversal is incorrect",1,result);
        Assert.assertEquals("Incorrect number of records in output",1,output.getRecords().size());
        Assert.assertEquals("Output record is incorrect",cleanedRead,output.getRecords().get(0));
    }

    @Test
    public void testMixedCleanedAndUncleanedNonOverlapping() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord[] sourceReads = new SAMRecord[] { ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read2",1,2,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read3",1,3,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read4",1,4,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read5",1,5,5) };

        SAMRecord cleanedRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,1);
        cleanedRead.setBaseQualities(getMockBaseQualityString((byte)1,cleanedRead.getReadLength()));
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader(cleanedRead);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:1-10", cleanedReads, output );

        int result = runWalkerOverReads(walker,sourceReads);

        Assert.assertEquals("Result of traversal is incorrect",1,result);
        Assert.assertEquals("Incorrect number of records in output",5,output.getRecords().size());
        for( int i = 0; i < sourceReads.length; i++ ) {
            if( i == 0 )
                Assert.assertEquals("Output record is incorrect",cleanedRead,output.getRecords().get(i));
            else
                Assert.assertEquals("Output record is incorrect",sourceReads[i],output.getRecords().get(i));
        }
    }

    @Test
    public void testAlternateReadOrder() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord[] sourceReads = new SAMRecord[] { ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read2",1,2,5),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read3",1,3,5) };

        SAMRecord cleanedRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,3,1);
        cleanedRead.setBaseQualities(getMockBaseQualityString((byte)1,cleanedRead.getReadLength()));
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader(cleanedRead);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:1-10", cleanedReads, output );

        int result = runWalkerOverReads(walker,sourceReads);

        Assert.assertEquals("Result of traversal is incorrect",1,result);
        Assert.assertEquals("Incorrect number of records in output",3,output.getRecords().size());
        Assert.assertEquals("Incorrect read at position 1",cleanedRead,output.getRecords().get(0));
        Assert.assertEquals("Incorrect read at position 2",sourceReads[1],output.getRecords().get(1));
        Assert.assertEquals("Incorrect read at position 3",sourceReads[2],output.getRecords().get(2));
    }

    @Test
    public void testReadOutsideInterval() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord sourceRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);

        SAMRecord cleanedRead = ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,5);
        cleanedRead.setBaseQualities(getMockBaseQualityString((byte)1,cleanedRead.getReadLength()));
        ArtificialSAMFileReader cleanedReads = new ArtificialSAMFileReader(cleanedRead);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:20-50", cleanedReads, output );

        int result = runWalkerOverReads( walker, sourceRead );

        Assert.assertEquals("Result of traversal is incorrect",0,result);
        Assert.assertEquals("Incorrect number of records in output",1,output.getRecords().size());
        Assert.assertEquals("Output record is incorrect",sourceRead,output.getRecords().get(0));
    }

    @Test
    public void testMultipleIntervals() {
        SAMFileHeader header = getMockSAMFileHeader();
        SAMRecord[] sourceReads = new SAMRecord[] { ArtificialSAMUtils.createArtificialRead(header,"read1",1,1,10),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read2",1,11,10),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read3",1,21,10),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read4",1,31,10),
                                                    ArtificialSAMUtils.createArtificialRead(header,"read5",1,41,10) };
        SAMRecord[] cleanedReads = new SAMRecord[sourceReads.length];
        for( int i = 0; i < sourceReads.length; i++ ) {
            try {
                cleanedReads[i] = (SAMRecord)sourceReads[i].clone();
                cleanedReads[i].setBaseQualities(getMockBaseQualityString((byte)1,cleanedReads[i].getReadLength()));
            }
            catch( CloneNotSupportedException ex ) {
                throw new StingException("Unable to clone samrecord", ex);
            }
        }

        ArtificialSAMFileReader cleanedReader = new ArtificialSAMFileReader(cleanedReads);

        ArtificialSAMFileWriter output = new ArtificialSAMFileWriter();
        CleanedReadInjector walker = createWalker( "chr1:11-20;chr1:31-40", cleanedReader, output );
        int result = runWalkerOverReads( walker, sourceReads );

        Assert.assertEquals("Result of traversal is incorrect",2,result);
        Assert.assertEquals("Incorrect number of records in output",5,output.getRecords().size());
        Assert.assertEquals("Incorrect read at position 1",sourceReads[0],output.getRecords().get(0));
        Assert.assertEquals("Incorrect read at position 2",cleanedReads[1],output.getRecords().get(1));
        Assert.assertEquals("Incorrect read at position 3",sourceReads[2],output.getRecords().get(2));
        Assert.assertEquals("Incorrect read at position 4",cleanedReads[3],output.getRecords().get(3));
        Assert.assertEquals("Incorrect read at position 5",sourceReads[4],output.getRecords().get(4));
    }

    private CleanedReadInjector createWalker( String intervals, ArtificialSAMFileReader cleanedReads, ArtificialSAMFileWriter output ) {
        CleanedReadInjector walker = new CleanedReadInjector();

        walker.intervalsSource = intervals;        
        walker.cleanedReadsSource = cleanedReads;
        walker.outputBAM = output;

        walker.initializeOutputStreams( new OutputTracker(null,null) );

        return walker;
    }

    private SAMFileHeader getMockSAMFileHeader() {
        return ArtificialSAMUtils.createArtificialSamHeader(2,0,247249719);
    }

    private Integer runWalkerOverReads( CleanedReadInjector walker, SAMRecord... reads ) {
        walker.initialize();
        Integer accum = walker.reduceInit();
        for( SAMRecord sourceRead: reads ) {
            Integer value = walker.map( null, sourceRead );
            accum = walker.reduce(value,accum);
        }
        walker.onTraversalDone(accum);
        return accum;
    }

    private byte[] getMockBaseQualityString( byte value, int length ) {
        byte[] baseQualities = new byte[length];
        Arrays.fill(baseQualities,value);
        return baseQualities;
    }

}
