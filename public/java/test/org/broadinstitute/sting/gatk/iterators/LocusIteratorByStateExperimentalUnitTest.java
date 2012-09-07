package org.broadinstitute.sting.gatk.iterators;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * testing of the experimental version of LocusIteratorByState
 */
public class LocusIteratorByStateExperimentalUnitTest extends BaseTest {
    private static SAMFileHeader header;
    private LocusIteratorByStateExperimental li;
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
    }

    private final LocusIteratorByStateExperimental makeLTBS(List<SAMRecord> reads, ReadProperties readAttributes) {
        return new LocusIteratorByStateExperimental(new FakeCloseableIterator<SAMRecord>(reads.iterator()), readAttributes, genomeLocParser, LocusIteratorByStateExperimental.sampleListForSAMWithoutReadGroups());
    }

    private static ReadProperties createTestReadProperties() {
        return createTestReadProperties(null);
    }

    private static ReadProperties createTestReadProperties( DownsamplingMethod downsamplingMethod ) {
        return new ReadProperties(
                Collections.<SAMReaderID>emptyList(),
                new SAMFileHeader(),
                false,
                SAMFileReader.ValidationStringency.STRICT,
                downsamplingMethod,
                new ValidationExclusion(),
                Collections.<ReadFilter>emptyList(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1
        );
    }

    private static class FakeCloseableIterator<T> implements CloseableIterator<T> {
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

    @Test
    public void testXandEQOperators() {
        final byte[] bases1 = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
        final byte[] bases2 = new byte[] {'A','A','A','C','A','A','A','A','A','C'};

        // create a test version of the Reads object
        ReadProperties readAttributes = createTestReadProperties();

        SAMRecord r1 = ArtificialSAMUtils.createArtificialRead(header,"r1",0,1,10);
        r1.setReadBases(bases1);
        r1.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r1.setCigarString("10M");

        SAMRecord r2 = ArtificialSAMUtils.createArtificialRead(header,"r2",0,1,10);
        r2.setReadBases(bases2);
        r2.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r2.setCigarString("3=1X5=1X");

        SAMRecord r3 = ArtificialSAMUtils.createArtificialRead(header,"r3",0,1,10);
        r3.setReadBases(bases2);
        r3.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
        r3.setCigarString("3=1X5M1X");

        SAMRecord r4  = ArtificialSAMUtils.createArtificialRead(header,"r4",0,1,10);
        r4.setReadBases(bases2);
        r4.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
        r4.setCigarString("10M");

        List<SAMRecord> reads = Arrays.asList(r1, r2, r3, r4);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads,readAttributes);

        while (li.hasNext()) {
            AlignmentContext context = li.next();
            ReadBackedPileup pileup = context.getBasePileup();
            Assert.assertEquals(pileup.depthOfCoverage(), 4);
        }
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

        List<SAMRecord> reads = Arrays.asList(before, during, after);

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

        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header, "indelOnly", 0, firstLocus, 76);
        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte) '@', 76));
        indelOnlyRead.setCigarString("76I");

        List<SAMRecord> reads = Arrays.asList(indelOnlyRead);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, readAttributes);

        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
        // and considers it to be an indel-containing read.
        Assert.assertTrue(li.hasNext(),"Should have found a whole-indel read in the normal base pileup without extended events enabled");
        AlignmentContext alignmentContext = li.next();
        Assert.assertEquals(alignmentContext.getLocation().getStart(), firstLocus, "Base pileup is at incorrect location.");
        ReadBackedPileup basePileup = alignmentContext.getBasePileup();
        Assert.assertEquals(basePileup.getReads().size(),1,"Pileup is of incorrect size");
        Assert.assertSame(basePileup.getReads().get(0), indelOnlyRead, "Read in pileup is incorrect");
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
        indelOnlyRead.setReadBases(Utils.dupBytes((byte) 'A', 76));
        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
        indelOnlyRead.setCigarString("76I");

        SAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',76));
        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
        fullMatchAfterIndel.setCigarString("75I1M");

        List<SAMRecord> reads = Arrays.asList(leadingRead, indelOnlyRead, fullMatchAfterIndel);

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

        Assert.assertEquals(numAlignmentContextsFound, 2, "Found incorrect number of alignment contexts");
    }

    /**
     * Test to make sure that reads supporting only an indel (example cigar string: 76I) are represented properly
     */
    @Test
    public void testWholeIndelReadRepresentedTest() {
        final int firstLocus = 44367788, secondLocus = firstLocus + 1;

        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,secondLocus,1);
        read1.setReadBases(Utils.dupBytes((byte) 'A', 1));
        read1.setBaseQualities(Utils.dupBytes((byte) '@', 1));
        read1.setCigarString("1I");

        List<SAMRecord> reads = Arrays.asList(read1);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getEventBases(), "A");
        }

        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,secondLocus,10);
        read2.setReadBases(Utils.dupBytes((byte) 'A', 10));
        read2.setBaseQualities(Utils.dupBytes((byte) '@', 10));
        read2.setCigarString("10I");

        reads = Arrays.asList(read2);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(reads, createTestReadProperties());

        while(li.hasNext()) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();
            Assert.assertTrue(pe.isBeforeInsertion());
            Assert.assertFalse(pe.isAfterInsertion());
            Assert.assertEquals(pe.getEventBases(), "AAAAAAAAAA");
        }
    }

    ////////////////////////////////////////////
    // comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////

    private static final int IS_BEFORE_DELETED_BASE_FLAG = 1;
    private static final int IS_BEFORE_DELETION_START_FLAG = 2;
    private static final int IS_AFTER_DELETED_BASE_FLAG = 4;
    private static final int IS_AFTER_DELETION_END_FLAG = 8;
    private static final int IS_BEFORE_INSERTION_FLAG = 16;
    private static final int IS_AFTER_INSERTION_FLAG = 32;
    private static final int IS_NEXT_TO_SOFTCLIP_FLAG = 64;

    private static class LIBSTest {


        final String cigar;
        final int readLength;
        final List<Integer> offsets;
        final List<Integer> flags;

        private LIBSTest(final String cigar, final int readLength, final List<Integer> offsets, final List<Integer> flags) {
            this.cigar = cigar;
            this.readLength = readLength;
            this.offsets = offsets;
            this.flags = flags;
        }
    }

    @DataProvider(name = "LIBSTest")
    public Object[][] createLIBSTestData() {
        return new Object[][]{
                {new LIBSTest("1I", 1, Arrays.asList(0), Arrays.asList(IS_BEFORE_INSERTION_FLAG))},
                {new LIBSTest("10I", 10, Arrays.asList(0), Arrays.asList(IS_BEFORE_INSERTION_FLAG))},
                {new LIBSTest("2M2I2M", 6, Arrays.asList(0,1,4,5), Arrays.asList(0,IS_BEFORE_INSERTION_FLAG,IS_AFTER_INSERTION_FLAG,0))},
                {new LIBSTest("2M2I", 4, Arrays.asList(0,1), Arrays.asList(0,IS_BEFORE_INSERTION_FLAG))},
                //TODO -- uncomment these when LIBS is fixed
                //{new LIBSTest("2I2M", 4, Arrays.asList(2,3), Arrays.asList(IS_AFTER_INSERTION_FLAG,0))},
                //{new LIBSTest("1I1M1D1M", 3, Arrays.asList(0,1), Arrays.asList(IS_AFTER_INSERTION_FLAG | IS_BEFORE_DELETION_START_FLAG | IS_BEFORE_DELETED_BASE_FLAG,IS_AFTER_DELETED_BASE_FLAG | IS_AFTER_DELETION_END_FLAG))},
                //{new LIBSTest("1S1I1M", 3, Arrays.asList(2), Arrays.asList(IS_AFTER_INSERTION_FLAG))},
                {new LIBSTest("1M2D2M", 3, Arrays.asList(0,1,2), Arrays.asList(IS_BEFORE_DELETION_START_FLAG | IS_BEFORE_DELETED_BASE_FLAG,IS_AFTER_DELETED_BASE_FLAG | IS_AFTER_DELETION_END_FLAG,0))},
                {new LIBSTest("1S1M", 2, Arrays.asList(1), Arrays.asList(IS_NEXT_TO_SOFTCLIP_FLAG))},
                {new LIBSTest("1M1S", 2, Arrays.asList(0), Arrays.asList(IS_NEXT_TO_SOFTCLIP_FLAG))},
                {new LIBSTest("1S1M1I", 3, Arrays.asList(1), Arrays.asList(IS_BEFORE_INSERTION_FLAG | IS_NEXT_TO_SOFTCLIP_FLAG))}
        };
    }

    @Test(dataProvider = "LIBSTest")
    public void testLIBS(LIBSTest params) {
        final int locus = 44367788;

        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, params.readLength);
        read.setReadBases(Utils.dupBytes((byte) 'A', params.readLength));
        read.setBaseQualities(Utils.dupBytes((byte) '@', params.readLength));
        read.setCigarString(params.cigar);

        // create the iterator by state with the fake reads and fake records
        li = makeLTBS(Arrays.asList(read), createTestReadProperties());

        int offset = 0;
        while ( li.hasNext() ) {
            AlignmentContext alignmentContext = li.next();
            ReadBackedPileup p = alignmentContext.getBasePileup();
            Assert.assertTrue(p.getNumberOfElements() == 1);
            PileupElement pe = p.iterator().next();

            final int flag = params.flags.get(offset);
            Assert.assertEquals(pe.isBeforeDeletedBase(), (flag & IS_BEFORE_DELETED_BASE_FLAG) != 0);
            Assert.assertEquals(pe.isBeforeDeletionStart(), (flag & IS_BEFORE_DELETION_START_FLAG) != 0);
            Assert.assertEquals(pe.isAfterDeletedBase(), (flag & IS_AFTER_DELETED_BASE_FLAG) != 0);
            Assert.assertEquals(pe.isAfterDeletionEnd(), (flag & IS_AFTER_DELETION_END_FLAG) != 0);
            Assert.assertEquals(pe.isBeforeInsertion(), (flag & IS_BEFORE_INSERTION_FLAG) != 0);
            Assert.assertEquals(pe.isAfterInsertion(), (flag & IS_AFTER_INSERTION_FLAG) != 0);
            Assert.assertEquals(pe.isNextToSoftClip(), (flag & IS_NEXT_TO_SOFTCLIP_FLAG) != 0);

            Assert.assertEquals(pe.getOffset(), params.offsets.get(offset).intValue());

            offset++;
        }
    }

    ////////////////////////////////////////////////
    // End comprehensive LIBS/PileupElement tests //
    ////////////////////////////////////////////////


    ///////////////////////////////////////
    // Read State Manager Tests          //
    ///////////////////////////////////////

    private class PerSampleReadStateManagerTest extends TestDataProvider {
        private List<Integer> readCountsPerAlignmentStart;
        private List<SAMRecord> reads;
        private List<ArrayList<LocusIteratorByStateExperimental.SAMRecordState>> recordStatesByAlignmentStart;
        private int removalInterval;

        public PerSampleReadStateManagerTest( List<Integer> readCountsPerAlignmentStart, int removalInterval ) {
            super(PerSampleReadStateManagerTest.class);

            this.readCountsPerAlignmentStart = readCountsPerAlignmentStart;
            this.removalInterval = removalInterval;

            reads = new ArrayList<SAMRecord>();
            recordStatesByAlignmentStart = new ArrayList<ArrayList<LocusIteratorByStateExperimental.SAMRecordState>>();

            setName(String.format("%s: readCountsPerAlignmentStart: %s  removalInterval: %d",
                                  getClass().getSimpleName(), readCountsPerAlignmentStart, removalInterval));
        }

        public void run() {
            LocusIteratorByStateExperimental libs = makeLTBS(new ArrayList<SAMRecord>(), createTestReadProperties());
            LocusIteratorByStateExperimental.ReadStateManager readStateManager =
                    libs.new ReadStateManager(new ArrayList<SAMRecord>().iterator());
            LocusIteratorByStateExperimental.ReadStateManager.PerSampleReadStateManager perSampleReadStateManager =
                    readStateManager.new PerSampleReadStateManager();

            makeReads();

            for ( ArrayList<LocusIteratorByStateExperimental.SAMRecordState> stackRecordStates : recordStatesByAlignmentStart ) {
                perSampleReadStateManager.addStatesAtNextAlignmentStart(stackRecordStates);
            }

            // read state manager should have the right number of reads
            Assert.assertEquals(reads.size(), perSampleReadStateManager.size());

            Iterator<SAMRecord> originalReadsIterator = reads.iterator();
            Iterator<LocusIteratorByStateExperimental.SAMRecordState> recordStateIterator = perSampleReadStateManager.iterator();
            int recordStateCount = 0;
            int numReadStatesRemoved = 0;

            // Do a first-pass validation of the record state iteration by making sure we get back everything we
            // put in, in the same order, doing any requested removals of read states along the way
            while ( recordStateIterator.hasNext() ) {
                LocusIteratorByStateExperimental.SAMRecordState readState = recordStateIterator.next();
                recordStateCount++;
                SAMRecord readFromPerSampleReadStateManager = readState.getRead();

                Assert.assertTrue(originalReadsIterator.hasNext());
                SAMRecord originalRead = originalReadsIterator.next();

                // The read we get back should be literally the same read in memory as we put in
                Assert.assertTrue(originalRead == readFromPerSampleReadStateManager);

                // If requested, remove a read state every removalInterval states
                if ( removalInterval > 0 && recordStateCount % removalInterval == 0 ) {
                    recordStateIterator.remove();
                    numReadStatesRemoved++;
                }
            }

            Assert.assertFalse(originalReadsIterator.hasNext());

            // If we removed any read states, do a second pass through the read states to make sure the right
            // states were removed
            if ( numReadStatesRemoved > 0 ) {
                Assert.assertEquals(perSampleReadStateManager.size(), reads.size() - numReadStatesRemoved);

                originalReadsIterator = reads.iterator();
                recordStateIterator = perSampleReadStateManager.iterator();
                int readCount = 0;
                int readStateCount = 0;

                // Match record states with the reads that should remain after removal
                while ( recordStateIterator.hasNext() ) {
                    LocusIteratorByStateExperimental.SAMRecordState readState = recordStateIterator.next();
                    readStateCount++;
                    SAMRecord readFromPerSampleReadStateManager = readState.getRead();

                    Assert.assertTrue(originalReadsIterator.hasNext());

                    SAMRecord originalRead = originalReadsIterator.next();
                    readCount++;

                    if ( readCount % removalInterval == 0 ) {
                        originalRead = originalReadsIterator.next(); // advance to next read, since the previous one should have been discarded
                        readCount++;
                    }

                    // The read we get back should be literally the same read in memory as we put in (after accounting for removals)
                    Assert.assertTrue(originalRead == readFromPerSampleReadStateManager);
                }

                Assert.assertEquals(readStateCount, reads.size() - numReadStatesRemoved);
            }

            // Allow memory used by this test to be reclaimed
            readCountsPerAlignmentStart = null;
            reads = null;
            recordStatesByAlignmentStart = null;
        }

        private void makeReads() {
            int alignmentStart = 1;

            for ( int readsThisStack : readCountsPerAlignmentStart ) {
                ArrayList<SAMRecord> stackReads = new ArrayList<SAMRecord>(ArtificialSAMUtils.createStackOfIdenticalArtificialReads(readsThisStack, header, "foo", 0, alignmentStart, MathUtils.randomIntegerInRange(50, 100)));
                ArrayList<LocusIteratorByStateExperimental.SAMRecordState> stackRecordStates = new ArrayList<LocusIteratorByStateExperimental.SAMRecordState>();

                for ( SAMRecord read : stackReads ) {
                    stackRecordStates.add(new LocusIteratorByStateExperimental.SAMRecordState(read));
                }

                reads.addAll(stackReads);
                recordStatesByAlignmentStart.add(stackRecordStates);
            }
        }
    }

    @DataProvider(name = "PerSampleReadStateManagerTestDataProvider")
    public Object[][] createPerSampleReadStateManagerTests() {
        for ( List<Integer> thisTestReadStateCounts : Arrays.asList( Arrays.asList(1),
                                                                     Arrays.asList(2),
                                                                     Arrays.asList(10),
                                                                     Arrays.asList(1, 1),
                                                                     Arrays.asList(2, 2),
                                                                     Arrays.asList(10, 10),
                                                                     Arrays.asList(1, 10),
                                                                     Arrays.asList(10, 1),
                                                                     Arrays.asList(1, 1, 1),
                                                                     Arrays.asList(2, 2, 2),
                                                                     Arrays.asList(10, 10, 10),
                                                                     Arrays.asList(1, 1, 1, 1, 1, 1),
                                                                     Arrays.asList(10, 10, 10, 10, 10, 10),
                                                                     Arrays.asList(1, 2, 10, 1, 2, 10)
                                                                   ) ) {

            for ( int removalInterval : Arrays.asList(0, 2, 3) ) {
                new PerSampleReadStateManagerTest(thisTestReadStateCounts, removalInterval);
            }
        }

        return PerSampleReadStateManagerTest.getTests(PerSampleReadStateManagerTest.class);
    }

    @Test(dataProvider = "PerSampleReadStateManagerTestDataProvider")
    public void runPerSampleReadStateManagerTest( PerSampleReadStateManagerTest test ) {
        logger.warn("Running test: " + test);

        test.run();
    }
}
