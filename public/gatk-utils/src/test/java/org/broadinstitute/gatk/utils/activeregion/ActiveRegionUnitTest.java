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

package org.broadinstitute.gatk.utils.activeregion;


// the imports for unit testing.


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public class ActiveRegionUnitTest extends BaseTest {
    private final static boolean DEBUG = false;
    private GenomeLocParser genomeLocParser;
    private IndexedFastaSequenceFile seq;
    private String contig;
    private int contigLength;

    @BeforeClass
    public void init() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        genomeLocParser = new GenomeLocParser(seq);
        contig = "1";
        contigLength = genomeLocParser.getContigInfo(contig).getSequenceLength();
    }

    @DataProvider(name = "ActionRegionCreationTest")
    public Object[][] makePollingData() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int size : Arrays.asList(1, 10, 100, 1000) ) {
                for ( final int ext : Arrays.asList(0, 1, 10, 100) ) {
                    for ( final boolean isActive : Arrays.asList(true, false) ) {
                        for ( final boolean addStates : Arrays.asList(true, false) ) {
                            List<ActivityProfileState> states = null;
                            if ( addStates ) {
                                states = new LinkedList<ActivityProfileState>();
                                for ( int i = start; i < start + size; i++ ) {
                                    states.add(new ActivityProfileState(genomeLocParser.createGenomeLoc(contig, i + start), isActive ? 1.0 : 0.0));
                                }
                            }
                            final GenomeLoc loc = genomeLocParser.createGenomeLoc(contig, start, start + size - 1);
                            tests.add(new Object[]{loc, states, isActive, ext});
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "ActionRegionCreationTest")
    public void testCreatingActiveRegions(final GenomeLoc loc, final List<ActivityProfileState> supportingStates, final boolean isActive, final int extension) {
        final ActiveRegion region = new ActiveRegion(loc, supportingStates, isActive, genomeLocParser, extension);
        Assert.assertEquals(region.getLocation(), loc);
        Assert.assertEquals(region.getExtendedLoc().getStart(), Math.max(loc.getStart() - extension, 1));
        Assert.assertEquals(region.getExtendedLoc().getStop(), Math.min(loc.getStop() + extension, contigLength));
        Assert.assertEquals(region.getReadSpanLoc().getStart(), Math.max(loc.getStart() - extension, 1));
        Assert.assertEquals(region.getReadSpanLoc().getStop(), Math.min(loc.getStop() + extension, contigLength));
        Assert.assertEquals(region.isActive(), isActive);
        Assert.assertEquals(region.getExtension(), extension);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getSupportingStates(), supportingStates == null ? Collections.emptyList() : supportingStates);
        Assert.assertNotNull(region.toString());

        assertGoodReferenceGetter(region.getActiveRegionReference(seq), region.getExtendedLoc(), 0);
        assertGoodReferenceGetter(region.getActiveRegionReference(seq, 10), region.getExtendedLoc(), 10);
        assertGoodReferenceGetter(region.getFullReference(seq), region.getReadSpanLoc(), 0);
        assertGoodReferenceGetter(region.getFullReference(seq, 10), region.getReadSpanLoc(), 10);
    }

    private void assertGoodReferenceGetter(final byte[] actualBytes, final GenomeLoc span, final int padding) {
        final int expectedStart = Math.max(span.getStart() - padding, 1);
        final int expectedStop = Math.min(span.getStop() + padding, contigLength);
        final byte[] expectedBytes = seq.getSubsequenceAt(span.getContig(), expectedStart, expectedStop).getBases();
        Assert.assertEquals(actualBytes, expectedBytes);
    }

    @DataProvider(name = "ActiveRegionReads")
    public Object[][] makeActiveRegionReads() {
        List<Object[]> tests = new ArrayList<Object[]>();
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        for ( final int start : Arrays.asList(1, 10, 100, contigLength - 10, contigLength - 1) ) {
            for ( final int readStartOffset : Arrays.asList(-100, -10, 0, 10, 100) ) {
                for ( final int readSize : Arrays.asList(10, 100, 1000) ) {
                    final GenomeLoc loc = genomeLocParser.createGenomeLocOnContig(contig, start, start + 10);

                    final int readStart = Math.max(start + readStartOffset, 1);
                    final int readStop = Math.min(readStart + readSize, contigLength);
                    final int readLength = readStop - readStart + 1;
                    if ( readLength > 0 ) {
                        GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, readStart, readLength);
                        final GenomeLoc readLoc = genomeLocParser.createGenomeLoc(read);
                        if ( readLoc.overlapsP(loc) )
                            tests.add(new Object[]{loc, read});
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "ActiveRegionReads")
    public void testActiveRegionReads(final GenomeLoc loc, final GATKSAMRecord read) throws Exception {
        final GenomeLoc expectedSpan = loc.union(genomeLocParser.createGenomeLoc(read));

        final ActiveRegion region = new ActiveRegion(loc, null, true, genomeLocParser, 0);
        final ActiveRegion region2 = new ActiveRegion(loc, null, true, genomeLocParser, 0);
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalExceptReads(region2));

        region.add(read);
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalExceptReads(region2));

        region.clearReads();
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalExceptReads(region2));

        region.addAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalExceptReads(region2));

        region.removeAll(Collections.<GATKSAMRecord>emptySet());
        Assert.assertEquals(region.getReads(), Collections.singletonList(read));
        Assert.assertEquals(region.size(), 1);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), expectedSpan);
        Assert.assertTrue(region.equalExceptReads(region2));

        region.removeAll(Collections.singleton(read));
        Assert.assertEquals(region.getReads(), Collections.emptyList());
        Assert.assertEquals(region.size(), 0);
        Assert.assertEquals(region.getExtendedLoc(), loc);
        Assert.assertEquals(region.getReadSpanLoc(), loc);
        Assert.assertTrue(region.equalExceptReads(region2));

        final GATKSAMRecord read2 = (GATKSAMRecord)read.clone();
        read2.setReadName(read.getReadName() + ".clone");

        for ( final GATKSAMRecord readToKeep : Arrays.asList(read, read2)) {
            region.addAll(Arrays.asList(read, read2));
            final GATKSAMRecord readToDiscard = readToKeep == read ? read2 : read;
            region.removeAll(Collections.singleton(readToDiscard));
            Assert.assertEquals(region.getReads(), Arrays.asList(readToKeep));
            Assert.assertEquals(region.size(), 1);
            Assert.assertEquals(region.getExtendedLoc(), loc);
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure bad inputs are properly detected
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "BadReadsTest")
    public Object[][] makeBadReadsTest() {
        List<Object[]> tests = new ArrayList<Object[]>();
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        tests.add(new Object[]{
                ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialSAMUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        tests.add(new Object[]{
                ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 10, 10),
                ArtificialSAMUtils.createArtificialRead(header, "read2", 1, 9, 10)});
        tests.add(new Object[]{
                ArtificialSAMUtils.createArtificialRead(header, "read1", 1, 10, 10),
                ArtificialSAMUtils.createArtificialRead(header, "read2", 0, 9, 10)});
        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = !DEBUG, dataProvider = "BadReadsTest", expectedExceptions = IllegalArgumentException.class)
    public void testBadReads(final GATKSAMRecord read1, final GATKSAMRecord read2) {
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(read1);
        final ActiveRegion region = new ActiveRegion(loc, null, true, genomeLocParser, 0);
        region.add(read1);
        region.add(read2);
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure we can properly cut up an active region based on engine intervals
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "SplitActiveRegion")
    public Object[][] makeSplitActiveRegion() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final GenomeLoc whole_span = genomeLocParser.createGenomeLoc("20", 1, 500);
        final GenomeLoc gl_before = genomeLocParser.createGenomeLoc("20", 1, 9);
        final GenomeLoc gl_after = genomeLocParser.createGenomeLoc("20", 250, 500);
        final GenomeLoc gl_diff_contig = genomeLocParser.createGenomeLoc("19", 40, 50);

        final int regionStart = 10;
        final int regionStop = 100;
        final GenomeLoc region = genomeLocParser.createGenomeLoc("20", regionStart, regionStop);

        for ( final GenomeLoc noEffect : Arrays.asList(whole_span) )
            tests.add(new Object[]{
                    region,
                    Arrays.asList(noEffect),
                    Arrays.asList(region)});

        for ( final GenomeLoc noOverlap : Arrays.asList(gl_before, gl_after, gl_diff_contig) )
            tests.add(new Object[]{
                    region,
                    Arrays.asList(noOverlap),
                    Arrays.asList()});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 5, 50)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", regionStart, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 50, 200)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 50, regionStop))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 40, 50)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 20, 30), genomeLocParser.createGenomeLoc("20", 40, 50)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 20, 30), genomeLocParser.createGenomeLoc("20", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 1, 30), genomeLocParser.createGenomeLoc("20", 40, 50)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", regionStart, 30), genomeLocParser.createGenomeLoc("20", 40, 50))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 1, 30), genomeLocParser.createGenomeLoc("20", 70, 200)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", regionStart, 30), genomeLocParser.createGenomeLoc("20", 70, regionStop))});

        tests.add(new Object[]{region,
                Arrays.asList(genomeLocParser.createGenomeLoc("20", 1, 30), genomeLocParser.createGenomeLoc("20", 40, 50), genomeLocParser.createGenomeLoc("20", 70, 200)),
                Arrays.asList(genomeLocParser.createGenomeLoc("20", regionStart, 30), genomeLocParser.createGenomeLoc("20", 40, 50), genomeLocParser.createGenomeLoc("20", 70, regionStop))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "SplitActiveRegion")
    public void testSplitActiveRegion(final GenomeLoc regionLoc, final List<GenomeLoc> intervalLocs, final List<GenomeLoc> expectedRegionLocs) {
        for ( final boolean addSubstates : Arrays.asList(true, false) ) {
            final List<ActivityProfileState> states;
            if ( addSubstates ) {
                states = new LinkedList<ActivityProfileState>();
                for ( int i = 0; i < regionLoc.size(); i++ )
                    states.add(new ActivityProfileState(genomeLocParser.createGenomeLoc(regionLoc.getContig(), regionLoc.getStart() + i), 1.0));
            } else {
                states = null;
            }

            final ActiveRegion region = new ActiveRegion(regionLoc, states, true, genomeLocParser, 0);
            final GenomeLocSortedSet intervals = new GenomeLocSortedSet(genomeLocParser,  intervalLocs);
            final List<ActiveRegion> regions = region.splitAndTrimToIntervals(intervals);

            Assert.assertEquals(regions.size(), expectedRegionLocs.size(), "Wrong number of split locations");
            for ( int i = 0; i < expectedRegionLocs.size(); i++ ) {
                final GenomeLoc expected = expectedRegionLocs.get(i);
                final ActiveRegion actual = regions.get(i);
                Assert.assertEquals(actual.getLocation(), expected, "Bad region after split");
                Assert.assertEquals(actual.isActive(), region.isActive());
                Assert.assertEquals(actual.getExtension(), region.getExtension());
            }
        }
    }

    // -----------------------------------------------------------------------------------------------
    //
    // Make sure we can properly cut up an active region based on engine intervals
    //
    // -----------------------------------------------------------------------------------------------

    @DataProvider(name = "TrimActiveRegionData")
    public Object[][] makeTrimActiveRegionData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // fully enclosed within active region
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 15, 16),
                genomeLocParser.createGenomeLoc("20", 15, 16), 0});

        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 10, 15),
                genomeLocParser.createGenomeLoc("20", 10, 15), 0});

        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 15, 20),
                genomeLocParser.createGenomeLoc("20", 15, 20), 0});

        // needs extra padding on the right
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 15, 25),
                genomeLocParser.createGenomeLoc("20", 15, 20), 5});

        // needs extra padding on the left
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 5, 15),
                genomeLocParser.createGenomeLoc("20", 10, 15), 5});

        // needs extra padding on both
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 7, 21),
                genomeLocParser.createGenomeLoc("20", 10, 20), 3});
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 9, 23),
                genomeLocParser.createGenomeLoc("20", 10, 20), 3});

        // desired span captures everything, so we're returning everything.  Tests that extension is set correctly
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 10, 20), 10,
                genomeLocParser.createGenomeLoc("20", 1, 50),
                genomeLocParser.createGenomeLoc("20", 10, 20), 10});

        // At the start of the chromosome, potentially a bit weird
        tests.add(new Object[]{
                genomeLocParser.createGenomeLoc("20", 1, 10), 10,
                genomeLocParser.createGenomeLoc("20", 1, 50),
                genomeLocParser.createGenomeLoc("20", 1, 10), 10});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "TrimActiveRegionData")
    public void testTrimActiveRegion(final GenomeLoc regionLoc, final int extension, final GenomeLoc desiredSpan, final GenomeLoc expectedActiveRegion, final int expectedExtension) {
        final ActiveRegion region = new ActiveRegion(regionLoc, Collections.<ActivityProfileState>emptyList(), true, genomeLocParser, extension);
        final ActiveRegion trimmed = region.trim(desiredSpan);
        Assert.assertEquals(trimmed.getLocation(), expectedActiveRegion, "Incorrect region");
        Assert.assertEquals(trimmed.getExtension(), expectedExtension, "Incorrect region");
    }
}