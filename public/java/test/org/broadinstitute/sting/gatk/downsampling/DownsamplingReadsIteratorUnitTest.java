package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIteratorAdapter;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collection;

public class DownsamplingReadsIteratorUnitTest {

    @Test
    public void testDownsamplingIteratorWithPositionalDownsampling() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        Collection<SAMRecord> reads = new ArrayList<SAMRecord>();

        reads.addAll(createStackOfIdenticalReads(3000, header, "foo", 0, 1, 100));
        reads.addAll(createStackOfIdenticalReads(3000, header, "foo", 0, 50, 100));

        StingSAMIterator iter = new DownsamplingReadsIterator(StingSAMIteratorAdapter.adapt(reads.iterator()), new PositionalDownsampler<SAMRecord>(1000));

        Assert.assertTrue(iter.hasNext());
        SAMRecord previous = iter.next();
        int count = 1;

        while ( iter.hasNext() ) {
            SAMRecord current = iter.next();
            Assert.assertTrue(previous.getAlignmentStart() <= current.getAlignmentStart() || ! previous.getReferenceIndex().equals(current.getReferenceIndex()));
            count++;
            previous = current;
        }

        Assert.assertEquals(count, 1000);
    }

    @Test
    public void testDownsamplingIteratorNoEffectiveDownsampling() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        Collection<SAMRecord> reads = new ArrayList<SAMRecord>();

        reads.addAll(createStackOfIdenticalReads(300, header, "foo", 0, 1, 100));
        reads.addAll(createStackOfIdenticalReads(300, header, "foo", 0, 50, 100));

        StingSAMIterator iter = new DownsamplingReadsIterator(StingSAMIteratorAdapter.adapt(reads.iterator()), new PositionalDownsampler<SAMRecord>(1000));

        Assert.assertTrue(iter.hasNext());
        SAMRecord previous = iter.next();
        int count = 1;

        while ( iter.hasNext() ) {
            SAMRecord current = iter.next();
            Assert.assertTrue(previous.getAlignmentStart() <= current.getAlignmentStart() || ! previous.getReferenceIndex().equals(current.getReferenceIndex()));
            count++;
            previous = current;
        }

        Assert.assertEquals(count, 600);
    }

    private ArrayList<SAMRecord> createStackOfIdenticalReads( int stackSize, SAMFileHeader header, String name, int refIndex, int alignmentStart, int length ) {
        ArrayList<SAMRecord> stack = new ArrayList<SAMRecord>(stackSize);
        for ( int i = 1; i <= stackSize; i++ ) {
            stack.add(ArtificialSAMUtils.createArtificialRead(header, name, refIndex, alignmentStart, length));
        }
        return stack;
    }
}
