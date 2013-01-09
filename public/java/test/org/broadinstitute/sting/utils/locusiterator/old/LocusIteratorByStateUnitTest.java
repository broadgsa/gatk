//package org.broadinstitute.sting.utils.locusiterator.old;
//
//import net.sf.samtools.*;
//import org.broadinstitute.sting.gatk.ReadProperties;
//import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
//import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
//import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
//import org.broadinstitute.sting.utils.NGSPlatform;
//import org.broadinstitute.sting.utils.Utils;
//import org.broadinstitute.sting.utils.locusiterator.LIBS_position;
//import org.broadinstitute.sting.utils.locusiterator.LocusIteratorByStateBaseTest;
//import org.broadinstitute.sting.utils.locusiterator.old.LocusIteratorByState;
//import org.broadinstitute.sting.utils.pileup.PileupElement;
//import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
//import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
//import org.broadinstitute.sting.utils.sam.GATKSAMReadGroupRecord;
//import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
//import org.testng.Assert;
//import org.testng.annotations.DataProvider;
//import org.testng.annotations.Test;
//
//import java.util.*;
//
///**
// * testing of the new (non-legacy) version of LocusIteratorByState
// */
//public class LocusIteratorByStateUnitTest extends LocusIteratorByStateBaseTest {
//
//    // TODO -- REMOVE ME WHEN LIBS IS FIXED
//    // TODO -- CURRENT CODE DOESN'T CORRECTLY COMPUTE THINGS LIKE BEFORE DELETION, AFTER INSERTION, ETC
//    private final static boolean ALLOW_BROKEN_LIBS_STATE = true;
//
//    protected LocusIteratorByState li;
//
//    @Test
//    public void testXandEQOperators() {
//        final byte[] bases1 = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
//        final byte[] bases2 = new byte[] {'A','A','A','C','A','A','A','A','A','C'};
//
//        // create a test version of the Reads object
//        ReadProperties readAttributes = createTestReadProperties();
//
//        SAMRecord r1 = ArtificialSAMUtils.createArtificialRead(header,"r1",0,1,10);
//        r1.setReadBases(bases1);
//        r1.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
//        r1.setCigarString("10M");
//
//        SAMRecord r2 = ArtificialSAMUtils.createArtificialRead(header,"r2",0,1,10);
//        r2.setReadBases(bases2);
//        r2.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
//        r2.setCigarString("3=1X5=1X");
//
//        SAMRecord r3 = ArtificialSAMUtils.createArtificialRead(header,"r3",0,1,10);
//        r3.setReadBases(bases2);
//        r3.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
//        r3.setCigarString("3=1X5M1X");
//
//        SAMRecord r4  = ArtificialSAMUtils.createArtificialRead(header,"r4",0,1,10);
//        r4.setReadBases(bases2);
//        r4.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
//        r4.setCigarString("10M");
//
//        List<SAMRecord> reads = Arrays.asList(r1, r2, r3, r4);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads,readAttributes);
//
//        while (li.hasNext()) {
//            AlignmentContext context = li.next();
//            ReadBackedPileup pileup = context.getBasePileup();
//            Assert.assertEquals(pileup.depthOfCoverage(), 4);
//        }
//    }
//
//    @Test
//    public void testIndelsInRegularPileup() {
//        final byte[] bases = new byte[] {'A','A','A','A','A','A','A','A','A','A'};
//        final byte[] indelBases = new byte[] {'A','A','A','A','C','T','A','A','A','A','A','A'};
//
//        // create a test version of the Reads object
//        ReadProperties readAttributes = createTestReadProperties();
//
//        SAMRecord before = ArtificialSAMUtils.createArtificialRead(header,"before",0,1,10);
//        before.setReadBases(bases);
//        before.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
//        before.setCigarString("10M");
//
//        SAMRecord during = ArtificialSAMUtils.createArtificialRead(header,"during",0,2,10);
//        during.setReadBases(indelBases);
//        during.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20,20,20});
//        during.setCigarString("4M2I6M");
//
//        SAMRecord after  = ArtificialSAMUtils.createArtificialRead(header,"after",0,3,10);
//        after.setReadBases(bases);
//        after.setBaseQualities(new byte[] {20,20,20,20,20,20,20,20,20,20});
//        after.setCigarString("10M");
//
//        List<SAMRecord> reads = Arrays.asList(before, during, after);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads,readAttributes);
//
//        boolean foundIndel = false;
//        while (li.hasNext()) {
//            AlignmentContext context = li.next();
//            ReadBackedPileup pileup = context.getBasePileup().getBaseFilteredPileup(10);
//            for (PileupElement p : pileup) {
//                if (p.isBeforeInsertion()) {
//                    foundIndel = true;
//                    Assert.assertEquals(p.getLengthOfImmediatelyFollowingIndel(), 2, "Wrong event length");
//                    Assert.assertEquals(p.getBasesOfImmediatelyFollowingInsertion(), "CT", "Inserted bases are incorrect");
//                    break;
//               }
//            }
//
//         }
//
//         Assert.assertTrue(foundIndel,"Indel in pileup not found");
//    }
//
//    @Test
//    public void testWholeIndelReadInIsolation() {
//        final int firstLocus = 44367789;
//
//        // create a test version of the Reads object
//        ReadProperties readAttributes = createTestReadProperties();
//
//        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header, "indelOnly", 0, firstLocus, 76);
//        indelOnlyRead.setReadBases(Utils.dupBytes((byte)'A',76));
//        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte) '@', 76));
//        indelOnlyRead.setCigarString("76I");
//
//        List<SAMRecord> reads = Arrays.asList(indelOnlyRead);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads, readAttributes);
//
//        // Traditionally, reads that end with indels bleed into the pileup at the following locus.  Verify that the next pileup contains this read
//        // and considers it to be an indel-containing read.
//        Assert.assertTrue(li.hasNext(),"Should have found a whole-indel read in the normal base pileup without extended events enabled");
//        AlignmentContext alignmentContext = li.next();
//        Assert.assertEquals(alignmentContext.getLocation().getStart(), firstLocus, "Base pileup is at incorrect location.");
//        ReadBackedPileup basePileup = alignmentContext.getBasePileup();
//        Assert.assertEquals(basePileup.getReads().size(),1,"Pileup is of incorrect size");
//        Assert.assertSame(basePileup.getReads().get(0), indelOnlyRead, "Read in pileup is incorrect");
//    }
//
//    /**
//     * Test to make sure that reads supporting only an indel (example cigar string: 76I) do
//     * not negatively influence the ordering of the pileup.
//     */
//    @Test
//    public void testWholeIndelRead() {
//        final int firstLocus = 44367788, secondLocus = firstLocus + 1;
//
//        SAMRecord leadingRead = ArtificialSAMUtils.createArtificialRead(header,"leading",0,firstLocus,76);
//        leadingRead.setReadBases(Utils.dupBytes((byte)'A',76));
//        leadingRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
//        leadingRead.setCigarString("1M75I");
//
//        SAMRecord indelOnlyRead = ArtificialSAMUtils.createArtificialRead(header,"indelOnly",0,secondLocus,76);
//        indelOnlyRead.setReadBases(Utils.dupBytes((byte) 'A', 76));
//        indelOnlyRead.setBaseQualities(Utils.dupBytes((byte)'@',76));
//        indelOnlyRead.setCigarString("76I");
//
//        SAMRecord fullMatchAfterIndel = ArtificialSAMUtils.createArtificialRead(header,"fullMatch",0,secondLocus,76);
//        fullMatchAfterIndel.setReadBases(Utils.dupBytes((byte)'A',76));
//        fullMatchAfterIndel.setBaseQualities(Utils.dupBytes((byte)'@',76));
//        fullMatchAfterIndel.setCigarString("75I1M");
//
//        List<SAMRecord> reads = Arrays.asList(leadingRead, indelOnlyRead, fullMatchAfterIndel);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads, createTestReadProperties());
//        int currentLocus = firstLocus;
//        int numAlignmentContextsFound = 0;
//
//        while(li.hasNext()) {
//            AlignmentContext alignmentContext = li.next();
//            Assert.assertEquals(alignmentContext.getLocation().getStart(),currentLocus,"Current locus returned by alignment context is incorrect");
//
//            if(currentLocus == firstLocus) {
//                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
//                Assert.assertEquals(readsAtLocus.size(),1,"Wrong number of reads at locus " + currentLocus);
//                Assert.assertSame(readsAtLocus.get(0),leadingRead,"leadingRead absent from pileup at locus " + currentLocus);
//            }
//            else if(currentLocus == secondLocus) {
//                List<GATKSAMRecord> readsAtLocus = alignmentContext.getBasePileup().getReads();
//                Assert.assertEquals(readsAtLocus.size(),2,"Wrong number of reads at locus " + currentLocus);
//                Assert.assertSame(readsAtLocus.get(0),indelOnlyRead,"indelOnlyRead absent from pileup at locus " + currentLocus);
//                Assert.assertSame(readsAtLocus.get(1),fullMatchAfterIndel,"fullMatchAfterIndel absent from pileup at locus " + currentLocus);
//            }
//
//            currentLocus++;
//            numAlignmentContextsFound++;
//        }
//
//        Assert.assertEquals(numAlignmentContextsFound, 2, "Found incorrect number of alignment contexts");
//    }
//
//    /**
//     * Test to make sure that reads supporting only an indel (example cigar string: 76I) are represented properly
//     */
//    @Test
//    public void testWholeIndelReadRepresentedTest() {
//        final int firstLocus = 44367788, secondLocus = firstLocus + 1;
//
//        SAMRecord read1 = ArtificialSAMUtils.createArtificialRead(header,"read1",0,secondLocus,1);
//        read1.setReadBases(Utils.dupBytes((byte) 'A', 1));
//        read1.setBaseQualities(Utils.dupBytes((byte) '@', 1));
//        read1.setCigarString("1I");
//
//        List<SAMRecord> reads = Arrays.asList(read1);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads, createTestReadProperties());
//
//        while(li.hasNext()) {
//            AlignmentContext alignmentContext = li.next();
//            ReadBackedPileup p = alignmentContext.getBasePileup();
//            Assert.assertTrue(p.getNumberOfElements() == 1);
//            PileupElement pe = p.iterator().next();
//            Assert.assertTrue(pe.isBeforeInsertion());
//            Assert.assertFalse(pe.isAfterInsertion());
//            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "A");
//        }
//
//        SAMRecord read2 = ArtificialSAMUtils.createArtificialRead(header,"read2",0,secondLocus,10);
//        read2.setReadBases(Utils.dupBytes((byte) 'A', 10));
//        read2.setBaseQualities(Utils.dupBytes((byte) '@', 10));
//        read2.setCigarString("10I");
//
//        reads = Arrays.asList(read2);
//
//        // create the iterator by state with the fake reads and fake records
//        li = makeLTBS(reads, createTestReadProperties());
//
//        while(li.hasNext()) {
//            AlignmentContext alignmentContext = li.next();
//            ReadBackedPileup p = alignmentContext.getBasePileup();
//            Assert.assertTrue(p.getNumberOfElements() == 1);
//            PileupElement pe = p.iterator().next();
//            Assert.assertTrue(pe.isBeforeInsertion());
//            Assert.assertFalse(pe.isAfterInsertion());
//            Assert.assertEquals(pe.getBasesOfImmediatelyFollowingInsertion(), "AAAAAAAAAA");
//        }
//    }
//
//    ////////////////////////////////////////////
//    // comprehensive LIBS/PileupElement tests //
//    ////////////////////////////////////////////
//
//    @DataProvider(name = "LIBSTest")
//    public Object[][] makeLIBSTest() {
//        final List<Object[]> tests = new LinkedList<Object[]>();
//
//        tests.add(new Object[]{new LIBSTest("1I", 1)});
//        tests.add(new Object[]{new LIBSTest("10I", 10)});
//        tests.add(new Object[]{new LIBSTest("2M2I2M", 6)});
//        tests.add(new Object[]{new LIBSTest("2M2I", 4)});
//        //TODO -- uncomment these when LIBS is fixed
//        //{new LIBSTest("2I2M", 4, Arrays.asList(2,3), Arrays.asList(IS_AFTER_INSERTION_FLAG,0))},
//        //{new LIBSTest("1I1M1D1M", 3, Arrays.asList(0,1), Arrays.asList(IS_AFTER_INSERTION_FLAG | IS_BEFORE_DELETION_START_FLAG | IS_BEFORE_DELETED_BASE_FLAG,IS_AFTER_DELETED_BASE_FLAG | IS_AFTER_DELETION_END_FLAG))},
//        //{new LIBSTest("1S1I1M", 3, Arrays.asList(2), Arrays.asList(IS_AFTER_INSERTION_FLAG))},
//        //{new LIBSTest("1M2D2M", 3)},
//        tests.add(new Object[]{new LIBSTest("1S1M", 2)});
//        tests.add(new Object[]{new LIBSTest("1M1S", 2)});
//        tests.add(new Object[]{new LIBSTest("1S1M1I", 3)});
//
//        return tests.toArray(new Object[][]{});
//
//        // TODO -- enable combinatorial tests here when LIBS is fixed
////        return createLIBSTests(
////                Arrays.asList(1, 10),
////                Arrays.asList(1, 2, 3));
//    }
//
//    @Test(dataProvider = "LIBSTest")
//    public void testLIBS(LIBSTest params) {
//        if ( params.getElements() == null || params.getElements().get(0).getOperator() == CigarOperator.I )
//            // TODO -- ENABLE ME WHEN LIBS IS FIXED
//            return;
//
//        // create the iterator by state with the fake reads and fake records
//        final GATKSAMRecord read = params.makeRead();
//        li = makeLTBS(Arrays.asList((SAMRecord)read), createTestReadProperties());
//        final LIBS_position tester = new LIBS_position(read);
//
//        int bpVisited = 0;
//        while ( li.hasNext() ) {
//            bpVisited++;
//
//            AlignmentContext alignmentContext = li.next();
//            ReadBackedPileup p = alignmentContext.getBasePileup();
//            Assert.assertTrue(p.getNumberOfElements() == 1);
//            PileupElement pe = p.iterator().next();
//
//            tester.stepForwardOnGenome();
//
//            if ( ! ALLOW_BROKEN_LIBS_STATE ) {
//                Assert.assertEquals(pe.isBeforeDeletedBase(), tester.isBeforeDeletedBase);
//                Assert.assertEquals(pe.isBeforeDeletionStart(), tester.isBeforeDeletionStart);
//                Assert.assertEquals(pe.isAfterDeletedBase(), tester.isAfterDeletedBase);
//                Assert.assertEquals(pe.isAfterDeletionEnd(), tester.isAfterDeletionEnd);
//                Assert.assertEquals(pe.isBeforeInsertion(), tester.isBeforeInsertion);
//                Assert.assertEquals(pe.isAfterInsertion(), tester.isAfterInsertion);
//                Assert.assertEquals(pe.isNextToSoftClip(), tester.isNextToSoftClip);
//            }
//
//            Assert.assertEquals(pe.getOffset(), tester.getCurrentReadOffset());
//        }
//
//        // min is one because always visit something, even for 10I reads
//        final int expectedBpToVisit = Math.max(read.getAlignmentEnd() - read.getAlignmentStart() + 1, 1);
//        Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
//    }
//
//    // ------------------------------------------------------------
//    //
//    // Tests for keeping reads
//    //
//    // ------------------------------------------------------------
//
//    @DataProvider(name = "LIBSKeepSubmittedReads")
//    public Object[][] makeLIBSKeepSubmittedReads() {
//        final List<Object[]> tests = new LinkedList<Object[]>();
//
//        for ( final boolean doSampling : Arrays.asList(true, false) ) {
//            for ( final int nReadsPerLocus : Arrays.asList(1, 10) ) {
//                for ( final int nLoci : Arrays.asList(1, 10, 25) ) {
//                    for ( final int nSamples : Arrays.asList(1, 2, 10) ) {
//                        for ( final boolean keepReads : Arrays.asList(true, false) ) {
//                            for ( final boolean grabReadsAfterEachCycle : Arrays.asList(true, false) ) {
////        for ( final int nReadsPerLocus : Arrays.asList(1) ) {
////            for ( final int nLoci : Arrays.asList(10) ) {
////                for ( final int nSamples : Arrays.asList(1) ) {
////                    for ( final boolean keepReads : Arrays.asList(true) ) {
////                        for ( final boolean grabReadsAfterEachCycle : Arrays.asList(true) ) {
//                                tests.add(new Object[]{nReadsPerLocus, nLoci, nSamples, keepReads, grabReadsAfterEachCycle, doSampling});
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        return tests.toArray(new Object[][]{});
//    }
//
//    @Test(enabled = true, dataProvider = "LIBSKeepSubmittedReads")
//    public void testLIBSKeepSubmittedReads(final int nReadsPerLocus,
//                                           final int nLoci,
//                                           final int nSamples,
//                                           final boolean keepReads,
//                                           final boolean grabReadsAfterEachCycle,
//                                           final boolean downsample) {
//        logger.warn(String.format("testLIBSKeepSubmittedReads %d %d %d %b %b %b", nReadsPerLocus, nLoci, nSamples, keepReads, grabReadsAfterEachCycle, downsample));
//        final int readLength = 10;
//
//        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 100000);
//        final List<String> samples = new ArrayList<String>(nSamples);
//        for ( int i = 0; i < nSamples; i++ ) {
//            final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("rg" + i);
//            final String sample = "sample" + i;
//            samples.add(sample);
//            rg.setSample(sample);
//            rg.setPlatform(NGSPlatform.ILLUMINA.getDefaultPlatform());
//            header.addReadGroup(rg);
//        }
//
//        final int maxCoveragePerSampleAtLocus = nReadsPerLocus * readLength / 2;
//        final int maxDownsampledCoverage = Math.max(maxCoveragePerSampleAtLocus / 2, 1);
//        final DownsamplingMethod downsampler = downsample
//                ? new DownsamplingMethod(DownsampleType.BY_SAMPLE, maxDownsampledCoverage, null, false)
//                : new DownsamplingMethod(DownsampleType.NONE, null, null, false);
//        final List<SAMRecord> reads = ArtificialSAMUtils.createReadStream(nReadsPerLocus, nLoci, header, 1, readLength);
//        li = new LocusIteratorByState(new FakeCloseableIterator<SAMRecord>(reads.iterator()),
//                createTestReadProperties(downsampler, keepReads),
//                genomeLocParser,
//                samples);
//
//        final Set<SAMRecord> seenSoFar = new HashSet<SAMRecord>();
//        final Set<SAMRecord> keptReads = new HashSet<SAMRecord>();
//        int bpVisited = 0;
//        while ( li.hasNext() ) {
//            bpVisited++;
//            final AlignmentContext alignmentContext = li.next();
//            final ReadBackedPileup p = alignmentContext.getBasePileup();
//
//            if ( downsample ) {
//                // just not a safe test
//                //Assert.assertTrue(p.getNumberOfElements() <= maxDownsampledCoverage * nSamples, "Too many reads at locus after downsampling");
//            } else {
//                final int minPileupSize = nReadsPerLocus * nSamples;
//                Assert.assertTrue(p.getNumberOfElements() >= minPileupSize);
//            }
//
//            seenSoFar.addAll(p.getReads());
//            if ( keepReads && grabReadsAfterEachCycle ) {
//                final List<SAMRecord> locusReads = li.transferReadsFromAllPreviousPileups();
//
//                // the number of reads starting here
//                int nReadsStartingHere = 0;
//                for ( final SAMRecord read : p.getReads() )
//                    if ( read.getAlignmentStart() == alignmentContext.getPosition() )
//                        nReadsStartingHere++;
//
//                if ( downsample )
//                    // with downsampling we might have some reads here that were downsampled away
//                    // in the pileup
//                    Assert.assertTrue(locusReads.size() >= nReadsStartingHere);
//                else
//                    Assert.assertEquals(locusReads.size(), nReadsStartingHere);
//                keptReads.addAll(locusReads);
//
//                // check that all reads we've seen so far are in our keptReads
//                for ( final SAMRecord read : seenSoFar ) {
//                    Assert.assertTrue(keptReads.contains(read), "A read that appeared in a pileup wasn't found in the kept reads: " + read);
//                }
//            }
//
//            if ( ! keepReads )
//                Assert.assertTrue(li.getReadsFromAllPreviousPileups().isEmpty(), "Not keeping reads but the underlying list of reads isn't empty");
//        }
//
//        if ( keepReads && ! grabReadsAfterEachCycle )
//            keptReads.addAll(li.transferReadsFromAllPreviousPileups());
//
//        if ( ! downsample ) { // downsampling may drop loci
//            final int expectedBpToVisit = nLoci + readLength - 1;
//            Assert.assertEquals(bpVisited, expectedBpToVisit, "Didn't visit the expected number of bp");
//        }
//
//        if ( keepReads ) {
//            // check we have the right number of reads
//            final int totalReads = nLoci * nReadsPerLocus * nSamples;
//            if ( ! downsample ) { // downsampling may drop reads
//                Assert.assertEquals(keptReads.size(), totalReads, "LIBS didn't keep the right number of reads during the traversal");
//
//                // check that the order of reads is the same as in our read list
//                for ( int i = 0; i < reads.size(); i++ ) {
//                    final SAMRecord inputRead = reads.get(i);
//                    final SAMRecord keptRead = reads.get(i);
//                    Assert.assertSame(keptRead, inputRead, "Input reads and kept reads differ at position " + i);
//                }
//            } else {
//                Assert.assertTrue(keptReads.size() <= totalReads, "LIBS didn't keep the right number of reads during the traversal");
//            }
//
//            // check uniqueness
//            final Set<String> readNames = new HashSet<String>();
//            for ( final SAMRecord read : keptReads ) {
//                Assert.assertFalse(readNames.contains(read.getReadName()), "Found duplicate reads in the kept reads");
//                readNames.add(read.getReadName());
//            }
//
//            // check that all reads we've seen are in our keptReads
//            for ( final SAMRecord read : seenSoFar ) {
//                Assert.assertTrue(keptReads.contains(read), "A read that appeared in a pileup wasn't found in the kept reads: " + read);
//            }
//        }
//    }
//}
