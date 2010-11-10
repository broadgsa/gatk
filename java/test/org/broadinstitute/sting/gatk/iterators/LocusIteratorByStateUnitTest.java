package org.broadinstitute.sting.gatk.iterators;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * testing of the LocusIteratorByState
 */
public class LocusIteratorByStateUnitTest extends BaseTest {

    private final int MAX_READS = 10;
    private static SAMFileHeader header;
    private LocusIteratorByState li;

    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    @Test
    public void testIndelBaseQualityFiltering() {
        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();
        JVMUtils.setFieldValue(JVMUtils.findField(ReadProperties.class,"generateExtendedEvents"),readAttributes,true);

        SAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
        before.setReadBases(bases);
        before.setBaseQualities(new byte[] {20,20,20,20,0,20,20,20,20,20});
        before.setCigarString("10M");

        SAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,2,10);
        during.setReadBases(bases);
        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20});
        during.setCigarString("4M1I6M");

        SAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,3,10);
        after.setReadBases(bases);
        after.setBaseQualities(new byte[] {20,20,0,20,20,20,20,20,20,20});
        after.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(before,during,after);

        // create the iterator by state with the fake reads and fake records
        li = new LocusIteratorByState(new FakeCloseableIterator<SAMRecord>(reads.iterator()),readAttributes,genomeLocParser);

        boolean foundExtendedEventPileup = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            if(!context.hasExtendedEventPileup())
                continue;

            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup().getBaseFilteredPileup(10);
            Assert.assertEquals(pileup.getLocation().getStart(), 5, "Extended event pileup at wrong location");
            Assert.assertEquals(pileup.size(), 3, "Pileup size is incorrect");

            foundExtendedEventPileup = true;
        }

        Assert.assertTrue(foundExtendedEventPileup,"Extended event pileup not found");
    }

    /**
     * Right now, the GATK's extended event pileup DOES NOT include reads which stop immediately before an insertion
     * but DOES include reads which stop immediately after an insertion.  This is almost certainly WRONG.  Eric is
     * figuring out the right way to handle this; in the meantime, adding this test to monitor that:
     *   A) the behavior is consistent
     *   B) so that we do end up with an automated test for this case when the model is fixed.
     */
    @Test
    public void testIndelPileupContainsAbuttingReads() {
        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] quals = new byte[] { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();
        JVMUtils.setFieldValue(JVMUtils.findField(ReadProperties.class,"generateExtendedEvents"),readAttributes,true);

        SAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
        before.setReadBases(bases);
        before.setBaseQualities(quals);
        before.setCigarString("10M");

        SAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,6,10);
        during.setReadBases(bases);
        during.setBaseQualities(quals);
        during.setCigarString("5M1I5M");

        SAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,11,10);
        after.setReadBases(bases);
        after.setBaseQualities(quals);
        after.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(before,during,after);

        // create the iterator by state with the fake reads and fake records
        li = new LocusIteratorByState(new FakeCloseableIterator<SAMRecord>(reads.iterator()),readAttributes,genomeLocParser);

        boolean foundExtendedEventPileup = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            if(!context.hasExtendedEventPileup())
                continue;

            Assert.assertEquals(context.getLocation().getStart(), 10, "Extended event pileup at wrong location");
            Assert.assertEquals(context.size(), 2, "Pileup size is incorrect");
            Assert.assertEquals(context.getExtendedEventPileup().getReads().get(0), during, "Read in pileup is incorrect");
            Assert.assertEquals(context.getExtendedEventPileup().getReads().get(1), after, "Read in pileup is incorrect");

            foundExtendedEventPileup = true;
        }

        Assert.assertTrue(foundExtendedEventPileup,"Extended event pileup not found");
    }

    private static ReadProperties createTestReadProperties() {
        return new ReadProperties(
                Collections.<SAMReaderID>emptyList(),
                new SAMFileHeader(),
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<SamRecordFilter>(),
                false,
                false
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
