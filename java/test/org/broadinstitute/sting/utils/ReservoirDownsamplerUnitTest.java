package org.broadinstitute.sting.utils;

import org.junit.Test;
import org.broadinstitute.sting.utils.sam.AlignmentStartComparator;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.*;

import junit.framework.Assert;

/**
 * Basic tests to prove the integrity of the reservoir downsampler.
 * At the moment, always run tests on SAM records as that's the task
 * for which the downsampler was conceived.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReservoirDownsamplerUnitTest {
    private static final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1,1,200);


    @Test
    public void testEmptyIterator() {
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(Collections.<SAMRecord>emptyList().iterator(),
                new AlignmentStartComparator(),1);
        Assert.assertFalse("Downsampler is not empty but should be.",downsampler.hasNext());
    }

    @Test
    public void testOneElementWithPoolSizeOne() {
        List<SAMRecord> reads = Collections.singletonList(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),1);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        Collection<SAMRecord> batchedReads = downsampler.next();
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertSame("Downsampler is returning an incorrect read",reads.get(0),batchedReads.iterator().next());
    }

    @Test
    public void testOneElementWithPoolSizeGreaterThanOne() {
        List<SAMRecord> reads = Collections.singletonList(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),5);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        Collection<SAMRecord> batchedReads = downsampler.next();
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertSame("Downsampler is returning an incorrect read",reads.get(0),batchedReads.iterator().next());

    }

    @Test
    public void testPoolFilledPartially() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),5);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",3,batchedReads.size());

        Assert.assertSame("Downsampler read 1 is incorrect",reads.get(0),batchedReads.get(0));
        Assert.assertSame("Downsampler read 2 is incorrect",reads.get(1),batchedReads.get(1));
        Assert.assertSame("Downsampler read 3 is incorrect",reads.get(2),batchedReads.get(2));
    }

    @Test
    public void testPoolFilledExactly() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),5);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",5,batchedReads.size());
        Assert.assertSame("Downsampler is returning an incorrect read",reads.get(0),batchedReads.iterator().next());

        Assert.assertSame("Downsampler read 1 is incorrect",reads.get(0),batchedReads.get(0));
        Assert.assertSame("Downsampler read 2 is incorrect",reads.get(1),batchedReads.get(1));
        Assert.assertSame("Downsampler read 3 is incorrect",reads.get(2),batchedReads.get(2));
        Assert.assertSame("Downsampler read 4 is incorrect",reads.get(3),batchedReads.get(3));
        Assert.assertSame("Downsampler read 5 is incorrect",reads.get(4),batchedReads.get(4));
    }

    @Test
    public void testLargerPileWithZeroElementPool() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),0);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",0,batchedReads.size());
        Assert.assertFalse("Downsampler is not empty but should be",downsampler.hasNext());
    }

    @Test
    public void testLargerPileWithSingleElementPool() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),1);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertTrue("Downsampler is returning a bad read.",reads.contains(batchedReads.get(0))) ;
        Assert.assertFalse("Downsampler is not empty but should be",downsampler.hasNext());
    }

    @Test
    public void testFillingAcrossLoci() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,2,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,2,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,3,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,3,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),5);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(0),batchedReads.get(0));

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",2,batchedReads.size());
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(1),batchedReads.get(0));
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(2),batchedReads.get(1));

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",2,batchedReads.size());
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(3),batchedReads.get(0));
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(4),batchedReads.get(1));

        Assert.assertFalse("Downsampler is not empty but should be",downsampler.hasNext());        
    }

    @Test
    public void testDownsamplingAcrossLoci() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,2,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,2,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,3,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,3,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(reads.iterator(),
                                                                                          new AlignmentStartComparator(),1);

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertEquals("Downsampler is returning an incorrect read.",reads.get(0),batchedReads.get(0));

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertTrue("Downsampler is returning an incorrect read.",batchedReads.get(0).equals(reads.get(1)) || batchedReads.get(0).equals(reads.get(2)));

        Assert.assertTrue("Downsampler is empty but shouldn't be",downsampler.hasNext());
        batchedReads = new ArrayList<SAMRecord>(downsampler.next());
        Assert.assertEquals("Downsampler is returning the wrong number of reads",1,batchedReads.size());
        Assert.assertTrue("Downsampler is returning an incorrect read.",batchedReads.get(0).equals(reads.get(3)) || batchedReads.get(0).equals(reads.get(4)));

        Assert.assertFalse("Downsampler is not empty but should be",downsampler.hasNext());
    }
}
