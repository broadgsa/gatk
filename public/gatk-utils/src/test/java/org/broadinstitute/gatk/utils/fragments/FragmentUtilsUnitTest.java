/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.fragments;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.gatk.utils.recalibration.EventType;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMReadGroupRecord;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Test routines for read-backed pileup.
 */
public class FragmentUtilsUnitTest extends BaseTest {
    private static SAMFileHeader header;
    private static GATKSAMReadGroupRecord rgForMerged;
    private final static boolean DEBUG = false;

    private class FragmentUtilsTest extends TestDataProvider {
        List<TestState> statesForPileup = new ArrayList<TestState>();
        List<TestState> statesForReads = new ArrayList<TestState>();

        private FragmentUtilsTest(String name, int readLen, int leftStart, int rightStart,
                                  boolean leftIsFirst, boolean leftIsNegative) {
            super(FragmentUtilsTest.class, String.format("%s-leftIsFirst:%b-leftIsNegative:%b", name, leftIsFirst, leftIsNegative));

            List<GATKSAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", readLen, leftStart, rightStart, leftIsFirst, leftIsNegative);
            GATKSAMRecord left = pair.get(0);
            GATKSAMRecord right = pair.get(1);

            for ( int pos = leftStart; pos < rightStart + readLen; pos++) {
                boolean posCoveredByLeft = pos >= left.getAlignmentStart() && pos <= left.getAlignmentEnd();
                boolean posCoveredByRight = pos >= right.getAlignmentStart() && pos <= right.getAlignmentEnd();

                if ( posCoveredByLeft || posCoveredByRight ) {
                    List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
                    List<Integer> offsets = new ArrayList<Integer>();

                    if ( posCoveredByLeft ) {
                        reads.add(left);
                        offsets.add(pos - left.getAlignmentStart());
                    }

                    if ( posCoveredByRight ) {
                        reads.add(right);
                        offsets.add(pos - right.getAlignmentStart());
                    }

                    boolean shouldBeFragment = posCoveredByLeft && posCoveredByRight;
                    ReadBackedPileup pileup = new ReadBackedPileupImpl(null, reads, offsets);
                    TestState testState = new TestState(shouldBeFragment ? 0 : 1, shouldBeFragment ? 1 : 0, pileup, null);
                    statesForPileup.add(testState);
                }

                TestState testState = left.getAlignmentEnd() >= right.getAlignmentStart() ? new TestState(0, 1, null, pair) : new TestState(2, 0, null, pair);
                statesForReads.add(testState);
            }
        }
    }

    private class TestState {
        int expectedSingletons, expectedPairs;
        ReadBackedPileup pileup;
        List<GATKSAMRecord> rawReads;

        private TestState(final int expectedSingletons, final int expectedPairs, final ReadBackedPileup pileup, final List<GATKSAMRecord> rawReads) {
            this.expectedSingletons = expectedSingletons;
            this.expectedPairs = expectedPairs;
            this.pileup = pileup;
            this.rawReads = rawReads;
        }
    }

    @DataProvider(name = "fragmentUtilsTest")
    public Object[][] createTests() {
        for ( boolean leftIsFirst : Arrays.asList(true, false) ) {
            for ( boolean leftIsNegative : Arrays.asList(true, false) ) {
                // Overlapping pair
                // ---->        [first]
                //   <---       [second]
                new FragmentUtilsTest("overlapping-pair", 10, 1, 5, leftIsFirst, leftIsNegative);

                // Non-overlapping pair
                // ---->
                //          <----
                new FragmentUtilsTest("nonoverlapping-pair", 10, 1, 15, leftIsFirst, leftIsNegative);
            }
        }

        return FragmentUtilsTest.getTests(FragmentUtilsTest.class);
    }

    @Test(enabled = !DEBUG, dataProvider = "fragmentUtilsTest")
    public void testAsPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            ReadBackedPileup rbp = testState.pileup;
            FragmentCollection<PileupElement> fp = FragmentUtils.create(rbp);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = !DEBUG, dataProvider = "fragmentUtilsTest")
    public void testAsListOfReadsFromPileup(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForPileup ) {
            FragmentCollection<GATKSAMRecord> fp = FragmentUtils.create(testState.pileup.getReads());
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = !DEBUG, dataProvider = "fragmentUtilsTest")
    public void testAsListOfReads(FragmentUtilsTest test) {
        for ( TestState testState : test.statesForReads ) {
            FragmentCollection<GATKSAMRecord> fp = FragmentUtils.create(testState.rawReads);
            Assert.assertEquals(fp.getOverlappingPairs().size(), testState.expectedPairs);
            Assert.assertEquals(fp.getSingletonReads().size(), testState.expectedSingletons);
        }
    }

    @Test(enabled = !DEBUG, expectedExceptions = IllegalArgumentException.class)
    public void testOutOfOrder() {
        final List<GATKSAMRecord> pair = ArtificialSAMUtils.createPair(header, "readpair", 100, 1, 50, true, true);
        final GATKSAMRecord left = pair.get(0);
        final GATKSAMRecord right = pair.get(1);
        final List<GATKSAMRecord> reads = Arrays.asList(right, left); // OUT OF ORDER!
        final List<Integer> offsets = Arrays.asList(0, 50);
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(null, reads, offsets);
        FragmentUtils.create(pileup); // should throw exception
    }

    @BeforeTest
    public void setup() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,1000);
        rgForMerged = new GATKSAMReadGroupRecord("RG1");
    }

    @DataProvider(name = "MergeFragmentsTest")
    public Object[][] createMergeFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) overlappingBaseQuals[i] = (byte)(i + 30);
            final GATKSAMRecord read1  = makeOverlappingRead(leftFlank, 20, overlappingBases, overlappingBaseQuals, "", 30, 1);
            final GATKSAMRecord read2  = makeOverlappingRead("", 20, overlappingBases, overlappingBaseQuals, rightFlank, 30, leftFlank.length() + 1);
            final GATKSAMRecord merged = makeOverlappingRead(leftFlank, 20, overlappingBases, overlappingBaseQuals, rightFlank, 30, 1);
            tests.add(new Object[]{"equalQuals", read1, read2, merged});

            // test that the merged read base quality is the
            tests.add(new Object[]{"lowQualLeft", modifyBaseQualities(read1, leftFlank.length(), overlapSize), read2, merged});
            tests.add(new Object[]{"lowQualRight", read1, modifyBaseQualities(read2, 0, overlapSize), merged});
        }

        return tests.toArray(new Object[][]{});
    }

    private GATKSAMRecord modifyBaseQualities(final GATKSAMRecord read, final int startOffset, final int length) throws Exception {
        final GATKSAMRecord readWithLowQuals = (GATKSAMRecord)read.clone();
        final byte[] withLowQuals = Arrays.copyOf(read.getBaseQualities(), read.getBaseQualities().length);
        for ( int i = startOffset; i < startOffset + length; i++ )
            withLowQuals[i] = (byte)(read.getBaseQualities()[i] + (i % 2 == 0 ? -1 : 0));
        readWithLowQuals.setBaseQualities(withLowQuals);
        return readWithLowQuals;
    }

    private GATKSAMRecord makeOverlappingRead(final String leftFlank, final int leftQual, final String overlapBases,
                                              final byte[] overlapQuals, final String rightFlank, final int rightQual,
                                              final int alignmentStart) {
        final String bases = leftFlank + overlapBases + rightFlank;
        final int readLength = bases.length();
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead", 0, alignmentStart, readLength);
        final byte[] leftQuals = Utils.dupBytes((byte) leftQual, leftFlank.length());
        final byte[] rightQuals = Utils.dupBytes((byte) rightQual, rightFlank.length());
        final byte[] quals = Utils.concat(leftQuals, overlapQuals, rightQuals);
        read.setCigarString(readLength + "M");
        read.setReadBases(bases.getBytes());
        for ( final EventType type : EventType.values() )
            read.setBaseQualities(quals, type);
        read.setReadGroup(rgForMerged);
        read.setMappingQuality(60);
        return read;
    }

    @Test(enabled = !DEBUG, dataProvider = "MergeFragmentsTest")
    public void testMergingTwoReads(final String name, final GATKSAMRecord read1, final GATKSAMRecord read2, final GATKSAMRecord expectedMerged) {
        final GATKSAMRecord actual = FragmentUtils.mergeOverlappingPairedFragments(read1, read2);

        if ( expectedMerged == null ) {
            Assert.assertNull(actual, "Expected reads not to merge, but got non-null result from merging");
        } else {
            Assert.assertTrue(actual.isStrandless(), "Merged reads should be strandless");
            Assert.assertNotNull(actual, "Expected reads to merge, but got null result from merging");
            // I really care about the bases, the quals, the CIGAR, and the read group tag
            Assert.assertEquals(actual.getCigarString(), expectedMerged.getCigarString());
            Assert.assertEquals(actual.getReadBases(), expectedMerged.getReadBases());
            Assert.assertEquals(actual.getReadGroup(), expectedMerged.getReadGroup());
            Assert.assertEquals(actual.getMappingQuality(), expectedMerged.getMappingQuality());
            for ( final EventType type : EventType.values() )
                Assert.assertEquals(actual.getBaseQualities(type), expectedMerged.getBaseQualities(type), "Failed base qualities for event type " + type);
        }
    }

    @Test(enabled = !DEBUG)
    public void testHardClippingBeforeMerge() {
        final String common = Utils.dupString("A", 10);
        final byte[] commonQuals = Utils.dupBytes((byte)30, common.length());
        final String adapter    = "NNNN";

        final GATKSAMRecord read1 = makeOverlappingRead(adapter, 30, common, commonQuals, "", 30, 10);
        final GATKSAMRecord read2 = makeOverlappingRead("", 30, common, commonQuals, adapter, 30, 10);
        final GATKSAMRecord expectedMerged = makeOverlappingRead("", 30, common, commonQuals, "", 30, 10);
        read1.setCigarString("4S" + common.length() + "M");
        read1.setProperPairFlag(true);
        read1.setReadPairedFlag(true);
        read1.setFirstOfPairFlag(true);
        read1.setReadNegativeStrandFlag(true);
        read1.setMateNegativeStrandFlag(false);
        read1.setMateAlignmentStart(read2.getAlignmentStart());
        read2.setCigarString(common.length() + "M4S");
        read2.setProperPairFlag(true);
        read2.setReadPairedFlag(true);
        read2.setFirstOfPairFlag(false);
        read2.setReadNegativeStrandFlag(false);
        read2.setMateNegativeStrandFlag(true);
        read2.setMateAlignmentStart(read1.getAlignmentStart());

        final int insertSize = common.length() - 1;
        read1.setInferredInsertSize(-insertSize);
        read2.setInferredInsertSize(insertSize);

        final GATKSAMRecord actual = FragmentUtils.mergeOverlappingPairedFragments(read1, read2);
        Assert.assertEquals(actual.getCigarString(), expectedMerged.getCigarString());
        Assert.assertEquals(actual.getReadBases(), expectedMerged.getReadBases());
        Assert.assertEquals(actual.getReadGroup(), expectedMerged.getReadGroup());
        Assert.assertEquals(actual.getMappingQuality(), expectedMerged.getMappingQuality());
        for ( final EventType type : EventType.values() )
            Assert.assertEquals(actual.getBaseQualities(type), expectedMerged.getBaseQualities(type), "Failed base qualities for event type " + type);
    }

    @Test(enabled = true)
    public void testHardClippingBeforeMergeResultingInCompletelyContainedSecondRead() {
        final String adapter    = "NNNN";

        final GATKSAMRecord read1 = makeOverlappingRead(adapter, 30, Utils.dupString("A", 10), Utils.dupBytes((byte)30, 10), "", 30, 10);
        final GATKSAMRecord read2 = makeOverlappingRead("", 30, Utils.dupString("A", 7), Utils.dupBytes((byte)30, 7), adapter, 30, 10);
        read1.setCigarString("4S10M");
        read1.setProperPairFlag(true);
        read1.setFirstOfPairFlag(true);
        read1.setReadNegativeStrandFlag(true);
        read1.setMateAlignmentStart(10);
        read2.setCigarString("7M4S");
        read2.setProperPairFlag(true);
        read2.setFirstOfPairFlag(false);
        read2.setReadNegativeStrandFlag(false);

        final int insertSize = 7 - 1;
        read1.setInferredInsertSize(insertSize);
        read2.setInferredInsertSize(-insertSize);

        final GATKSAMRecord actual = FragmentUtils.mergeOverlappingPairedFragments(read1, read2);
        Assert.assertNull(actual);
    }

    @DataProvider(name = "MergeFragmentsOffContig")
    public Object[][] makeMergeFragmentsOffContig() throws Exception {
        List<Object[]> tests = new ArrayList<>();

        for ( final int pre1 : Arrays.asList(0, 50)) {
            for ( final int post1 : Arrays.asList(0, 50)) {
                for ( final int pre2 : Arrays.asList(0, 50)) {
                    for ( final int post2 : Arrays.asList(0, 50)) {
                        tests.add(new Object[]{pre1, post1, pre2, post2});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MergeFragmentsOffContig")
    public void testMergeFragmentsOffContig(final int pre1, final int post1, final int pre2, final int post2) {
        final int contigSize = 10;
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 0, contigSize);

        final GATKSAMRecord read1 = createReadOffContig(header, false, pre1, post1);
        final GATKSAMRecord read2 = createReadOffContig(header, true, pre2, post2);

        final GATKSAMRecord merged = FragmentUtils.mergeOverlappingPairedFragments(read1, read2);
    }

    private GATKSAMRecord createReadOffContig(final SAMFileHeader header, final boolean negStrand, final int pre, final int post) {
        final int contigLen = header.getSequence(0).getSequenceLength();
        final int readLen = pre + contigLen + post;
        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 1, readLen);
        read.setAlignmentStart(1);
        read.setCigar(TextCigarCodec.decode(pre + "S" + contigLen + "M" + post + "S"));
        read.setBaseQualities(Utils.dupBytes((byte) 30, readLen));
        read.setReadBases(Utils.dupBytes((byte)'A', readLen));
        read.setMappingQuality(60);
        read.setMateAlignmentStart(1);
        read.setProperPairFlag(true);
        read.setReadPairedFlag(true);
        read.setInferredInsertSize(30);
        read.setReadNegativeStrandFlag(negStrand);
        read.setMateNegativeStrandFlag(! negStrand);
        read.setReadGroup(new GATKSAMReadGroupRecord("foo"));
        return read;
    }


    private static final byte highQuality = 30;
    private static final byte overlappingQuality = 20;

    @DataProvider(name = "AdjustFragmentsTest")
    public Object[][] createAdjustFragmentsTest() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        final String leftFlank = "CCC";
        final String rightFlank = "AAA";
        final String allOverlappingBases = "ACGTACGTGGAACCTTAG";
        for ( int overlapSize = 1; overlapSize < allOverlappingBases.length(); overlapSize++ ) {
            final String overlappingBases = allOverlappingBases.substring(0, overlapSize);
            final byte[] overlappingBaseQuals = new byte[overlapSize];
            for ( int i = 0; i < overlapSize; i++ ) overlappingBaseQuals[i] = highQuality;
            final GATKSAMRecord read1  = makeOverlappingRead(leftFlank, highQuality, overlappingBases, overlappingBaseQuals, "", highQuality, 1);
            final GATKSAMRecord read2  = makeOverlappingRead("", highQuality, overlappingBases, overlappingBaseQuals, rightFlank, highQuality, leftFlank.length() + 1);
            tests.add(new Object[]{read1, read2, overlapSize});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "AdjustFragmentsTest")
    public void testAdjustingTwoReads(final GATKSAMRecord read1, final GATKSAMRecord read2, final int overlapSize) {
        FragmentUtils.adjustQualsOfOverlappingPairedFragments(read1, read2);

        for ( int i = 0; i < read1.getReadLength() - overlapSize; i++ )
            Assert.assertEquals(read1.getBaseQualities()[i], highQuality);
        for ( int i = read1.getReadLength() - overlapSize; i < read1.getReadLength(); i++ )
            Assert.assertEquals(read1.getBaseQualities()[i], overlappingQuality);

        for ( int i = 0; i < overlapSize; i++ )
            Assert.assertEquals(read2.getBaseQualities()[i], overlappingQuality);
        for ( int i = overlapSize; i < read2.getReadLength(); i++ )
            Assert.assertEquals(read2.getBaseQualities()[i], highQuality);
    }
}
