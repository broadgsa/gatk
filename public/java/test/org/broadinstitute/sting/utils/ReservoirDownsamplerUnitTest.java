package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.*;

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
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(1);
        Assert.assertTrue(downsampler.isEmpty(),"Downsampler is not empty but should be.");
    }

    @Test
    public void testOneElementWithPoolSizeOne() {
        List<GATKSAMRecord> reads = Collections.singletonList(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(1);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        Collection<SAMRecord> batchedReads = downsampler.getDownsampledContents();
        Assert.assertEquals(batchedReads.size(), 1, "Downsampler is returning the wrong number of reads");
        Assert.assertSame(batchedReads.iterator().next(), reads.get(0), "Downsampler is returning an incorrect read");
    }

    @Test
    public void testOneElementWithPoolSizeGreaterThanOne() {
        List<GATKSAMRecord> reads = Collections.singletonList(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(5);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        Collection<SAMRecord> batchedReads = downsampler.getDownsampledContents();
        Assert.assertEquals(batchedReads.size(), 1, "Downsampler is returning the wrong number of reads");
        Assert.assertSame(batchedReads.iterator().next(), reads.get(0), "Downsampler is returning an incorrect read");

    }

    @Test
    public void testPoolFilledPartially() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(5);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 3, "Downsampler is returning the wrong number of reads");

        Assert.assertSame(batchedReads.get(0), reads.get(0), "Downsampler read 1 is incorrect");
        Assert.assertSame(batchedReads.get(1), reads.get(1), "Downsampler read 2 is incorrect");
        Assert.assertSame(batchedReads.get(2), reads.get(2), "Downsampler read 3 is incorrect");
    }

    @Test
    public void testPoolFilledExactly() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(5);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 5, "Downsampler is returning the wrong number of reads");
        Assert.assertSame(batchedReads.iterator().next(), reads.get(0), "Downsampler is returning an incorrect read");

        Assert.assertSame(batchedReads.get(0), reads.get(0), "Downsampler read 1 is incorrect");
        Assert.assertSame(batchedReads.get(1), reads.get(1), "Downsampler read 2 is incorrect");
        Assert.assertSame(batchedReads.get(2), reads.get(2), "Downsampler read 3 is incorrect");
        Assert.assertSame(batchedReads.get(3), reads.get(3), "Downsampler read 4 is incorrect");
        Assert.assertSame(batchedReads.get(4), reads.get(4), "Downsampler read 5 is incorrect");
    }

    @Test
    public void testLargerPileWithZeroElementPool() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(0);
        downsampler.addAll(reads);

        Assert.assertTrue(downsampler.isEmpty(),"Downsampler isn't empty but should be");
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 0, "Downsampler is returning the wrong number of reads");
    }

    @Test
    public void testLargerPileWithSingleElementPool() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,1,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(1);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 1, "Downsampler is returning the wrong number of reads");
        Assert.assertTrue(reads.contains(batchedReads.get(0)),"Downsampler is returning a bad read.");
    }

    @Test
    public void testFillingAcrossLoci() {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read1",0,1,76));
        ReservoirDownsampler<SAMRecord> downsampler = new ReservoirDownsampler<SAMRecord>(5);
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        List<SAMRecord> batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 1, "Downsampler is returning the wrong number of reads");
        Assert.assertEquals(batchedReads.get(0), reads.get(0), "Downsampler is returning an incorrect read.");

        reads.clear();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read2",0,2,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read3",0,2,76));

        downsampler.clear();
        downsampler.addAll(reads);

        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 2, "Downsampler is returning the wrong number of reads");
        Assert.assertEquals(batchedReads.get(0), reads.get(0), "Downsampler is returning an incorrect read.");
        Assert.assertEquals(batchedReads.get(1), reads.get(1), "Downsampler is returning an incorrect read.");

        reads.clear();
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read4",0,3,76));
        reads.add(ArtificialSAMUtils.createArtificialRead(header,"read5",0,3,76));

        downsampler.clear();
        downsampler.addAll(reads);
                
        Assert.assertFalse(downsampler.isEmpty(),"Downsampler is empty but shouldn't be");
        batchedReads = new ArrayList<SAMRecord>(downsampler.getDownsampledContents());
        Assert.assertEquals(batchedReads.size(), 2, "Downsampler is returning the wrong number of reads");
        Assert.assertEquals(batchedReads.get(0), reads.get(0), "Downsampler is returning an incorrect read.");
        Assert.assertEquals(batchedReads.get(1), reads.get(1), "Downsampler is returning an incorrect read.");
    }

}
