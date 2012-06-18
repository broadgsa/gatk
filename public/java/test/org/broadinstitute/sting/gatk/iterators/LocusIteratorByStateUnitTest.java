package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * testing of the LocusIteratorByState
 */
public class LocusIteratorByStateUnitTest extends BaseTest {
    private static SAMFileHeader header;
    private LocusIteratorByState li;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    private final LocusIteratorByState makeLTBS(List<SAMRecord> reads, ReadProperties readAttributes) {
        return new LocusIteratorByState(new FakeCloseableIterator<SAMRecord>(reads.iterator()), readAttributes, genomeLocParser, LocusIteratorByState.sampleListForSAMWithoutReadGroups());
    }

    @Test
    public void testIndelsInRegularPileup() {
        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] indelBases = new byte[] {'A','A','A','A','C','T','A','A','A','A','A','A'};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
        before.setReadBases(bases);
        before.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        before.setCigarString("10M");

        SAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,2,10);
        during.setReadBases(indelBases);
        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        during.setCigarString("4M2I6M");

        SAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,3,10);
        after.setReadBases(bases);
        after.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        after.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(before,during,after);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads,readAttributes);

        boolean foundIndel = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup().getBaseFilteredPileup(10);
            for (PileupElement p : pileup) {
                if (p.isBeforeInsertion()) {
                    foundIndel = true;
                    Assert.assertEquals(p.getEventLength(), 2, "Wrong event length");
                    Assert.assertEquals(p.getEventBases(), "CT", "Inserted bases are incorrect");
                    break;
               }
            }

         }

         Assert.assertTrue(foundIndel,"Indel in pileup not found");
    }

    @Test
    public void testWholeIndelReadInIsolation() {
        final int firstLocus = 44367789;

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,firstLocus,76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        List<SAMRecord> reads = Arrays.asList(indelOnlyRead);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, readAttributes);

        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
        // and considers it to be an indel-containing read.
        Assert.assertTrue(li.hasNext(),"Should have found a whole-indel read in the normal base pileup without extended events enabled");
        AlignmentContext alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(),firstLocus,"Base pileup is at incorrect location.");
        ReadBackedPileup basePileup = alignmentContext.getBasePileup();
        Assert.assertEquals(basePileup.getReads().size(),1,"Pileup is of incorrect size");
        Assert.assertSame(basePileup.getReads().get(0),indelOnlyRead,"Read in pileup is incorrect");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test
    public void testWholeIndelRead() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        SAMRecord leadingRead = ArtificialSAMUtils.createArtificialRead(header,"leading",0,firstLocus,76);
        leadingRead.setReadBases(Utils.dupBytes((byte)'A',76));
        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        leadingRead.setCigarString("1M75I");

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        SAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',76));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
        fullMatchAfterIndel.setCigarString("75I1M");

        List<SAMRecord> reads = Arrays.asList(leadingRead,indelOnlyRead,fullMatchAfterIndel);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());
        int currentLocus = firstLocus;
        int numAlignmentContextsFound = 0;

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            Assert.assertEquals(alignmentContext.getLocation().getStart(),currentLocus,"Current locus returned by alignment context is incorrect");

            if(currentLocus == firstLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + currentLocus);
            }
            else if(currentLocus == secondLocus) {
                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
                Assert.assertEquals(readsAtLocus.size(),2,"Wrong number of reads at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(0),indelOnlyRead,"indelOnlyRead absent from pileup at locus " + currentLocus);
                Assert.assertSame(readsAtLocus.get(1),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at locus " + currentLocus);
            }

            currentLocus++;
            numAlignmentContextsFound++;
        }

        Assert.assertEquals(numAlignmentContextsFound,2,"Found incorrect number of alignment contexts");
    }

    private static ReadProperties createTestReadProperties() {
        return new ReadProperties(
                Collections.<SAMReaderID>emptyList(),
                new SAMFileHeader(),
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                new ValidationExclusion(),
                Collections.<ReadFilter>emptyList(),
                false,
                BAQ.CalculationMode.OFF,
                BAQ.QualityMode.DONT_MODIFY,
                null, // no BAQ
                null, // no BQSR
                (byte) -1
        );
    }
}

class FakeCloseableIterator<T> implements CloseableIterator<T> {
    Iterator<T> iterator;

    public FakeCloseableIterator(Iterator<T> it) {
        iterator = it;
    }

    @Override
    public void close() {
        return;
    }

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public T next() {
        return iterator.next();
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Don't remove!");
    }
}
