package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.Utils;
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
        li = makeLTBS(reads,readAttributes);

        boolean foundExtendedEventPileup = false;
        while (li.hasNext()) {
            AlignmentContext context = li.next();
            if(!context.hasExtendedEventPileup())
                continue;

            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup().getBaseFilteredPileup(10);
            Assert.assertEquals(pileup.getLocation().getStart(), 5, "Extended event pileup at wrong location");
            Assert.assertEquals(pileup.getNumberOfElements(), 3, "Pileup size is incorrect");

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
        li = makeLTBS(reads,readAttributes);

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

        // Turn on extended events, and make sure the event is found.
        JVMUtils.setFieldValue(JVMUtils.findField(ReadProperties.class,"generateExtendedEvents"),readAttributes,true);
        li = makeLTBS(reads, readAttributes);

        Assert.assertTrue(li.hasNext(),"LocusIteratorByState with extended events should contain exactly one pileup");
        alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(),firstLocus-1,"Extended event pileup is at incorrect location.");
        ReadBackedExtendedEventPileup extendedEventPileup = alignmentContext.getExtendedEventPileup();
        Assert.assertEquals(extendedEventPileup.getReads().size(),1,"Pileup is of incorrect size");
        Assert.assertSame(extendedEventPileup.getReads().get(0),indelOnlyRead,"Read in pileup is incorrect");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test
    public void testWholeIndelReadWithoutExtendedEvents() {
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

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
     * not negatively influence the ordering of the pileup.
     */
    @Test
    public void testWholeIndelReadWithExtendedEvents() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();
        JVMUtils.setFieldValue(JVMUtils.findField(ReadProperties.class,"generateExtendedEvents"),readAttributes,true);

        SAMRecord leadingRead = ArtificialSAMUtils.createArtificialRead(header,"leading",0,firstLocus,76);
        leadingRead.setReadBases(Utils.dupBytes((byte)'A',76));
        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        leadingRead.setCigarString("1M75I");

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        SAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,1);
        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',1));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',1));
        fullMatchAfterIndel.setCigarString("1M");

        List<SAMRecord> reads = Arrays.asList(leadingRead,indelOnlyRead,fullMatchAfterIndel);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads,readAttributes);

        Assert.assertTrue(li.hasNext(),"Missing first locus at " + firstLocus);
        AlignmentContext alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(),firstLocus,"Incorrect locus at this position; should be " + firstLocus);
        List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
        Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + firstLocus);
        Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + firstLocus);

        Assert.assertTrue(li.hasNext(),"Missing extended event at " + firstLocus);
        alignmentContext = li.next();        
        Assert.assertEquals(alignmentContext.getLocation().getStart(),firstLocus,"Incorrect extended event locus at this position; should be " + firstLocus);
        readsAtLocus = alignmentContext.getExtendedEventPileup().getReads();
        Assert.assertEquals(readsAtLocus.size(),3,"Wrong number of reads at extended event locus " + firstLocus);
        Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at extended event locus " + firstLocus);
        Assert.assertSame(readsAtLocus.get(1),indelOnlyRead,"indelOnlyRead absent from pileup at extended event locus " + firstLocus);
        // Weird, but as above, reads immediately after the indel are included in the extended event pileup
        Assert.assertSame(readsAtLocus.get(2),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at extended event locus " + firstLocus);

        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
        // and considers it to be an indel-containing read.        
        Assert.assertTrue(li.hasNext(),"Missing base pileup at " + secondLocus);
        alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(),secondLocus,"Incorrect extended event locus at this position; should be " + secondLocus);
        readsAtLocus = alignmentContext.getBasePileup().getReads();
        Assert.assertEquals(readsAtLocus.size(),3,"Wrong number of reads at extended event locus " + secondLocus);
        Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at extended event locus " + secondLocus);
        Assert.assertSame(readsAtLocus.get(1),indelOnlyRead,"indelOnlyRead absent from pileup at extended event locus " + secondLocus);
        // Weird, but as above, reads immediately after the indel are included in the extended event pileup
        Assert.assertSame(readsAtLocus.get(2),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at extended event locus " + secondLocus);

        Assert.assertFalse(li.hasNext(),"Too many alignment contexts");
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
                Collections.<ReadFilter>emptyList(),
                false,
                false,
                BAQ.CalculationMode.OFF,
                BAQ.QualityMode.DONT_MODIFY,
                null, // no BAQ
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
