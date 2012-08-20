package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.List;

public class FractionalDownsamplerUnitTest {

    @Test
    public void test100PercentInclusion() {
        FractionalDownsampler<SAMRecord> downsampler = new FractionalDownsampler<SAMRecord>(1.0);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        downsampler.submit(createRandomReads(1000, header, "foo", 0, 100000, 500));
        downsampler.signalEndOfInput();

        List<SAMRecord> downsampledReads = downsampler.consumeDownsampledItems();

        Assert.assertTrue(downsampledReads.size() == 1000);
    }

    @Test
    public void test0PercentInclusion() {
        FractionalDownsampler<SAMRecord> downsampler = new FractionalDownsampler<SAMRecord>(0.0);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        downsampler.submit(createRandomReads(1000, header, "foo", 0, 100000, 500));
        downsampler.signalEndOfInput();

        List<SAMRecord> downsampledReads = downsampler.consumeDownsampledItems();

        Assert.assertTrue(downsampledReads.isEmpty());
    }

    @Test
    public void test50PercentInclusion() {
        FractionalDownsampler<SAMRecord> downsampler = new FractionalDownsampler<SAMRecord>(0.5);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000);

        downsampler.submit(createRandomReads(5000, header, "foo", 0, 100000, 500));
        downsampler.signalEndOfInput();

        List<SAMRecord> downsampledReads = downsampler.consumeDownsampledItems();

        Assert.assertTrue(downsampledReads.size() >= 2000 && downsampledReads.size() <= 3000);
    }

    private List<SAMRecord> createRandomReads( int numReads, SAMFileHeader header, String name, int contigIndex, int maxAlignmentStart, int maxLength ) {
        List<SAMRecord> reads = new ArrayList<SAMRecord>(numReads);

        for ( int i = 1; i <= numReads; i++ ) {
            reads.add(ArtificialSAMUtils.createArtificialRead(header, name, contigIndex,
                                                              GenomeAnalysisEngine.getRandomGenerator().nextInt(maxAlignmentStart) + 1,
                                                              GenomeAnalysisEngine.getRandomGenerator().nextInt(maxLength) + 1));
        }

        return reads;
    }
}
