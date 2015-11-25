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

package org.broadinstitute.gatk.engine.traversals;

import com.google.java.contract.PreconditionError;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.datasources.reads.*;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.broadinstitute.gatk.utils.sam.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.engine.executive.WindowMaker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 11/13/12
 * Time: 2:47 PM
 *
 * Test the Active Region Traversal Contract
 * http://iwww.broadinstitute.org/gsa/wiki/index.php/Active_Region_Traversal_Contract
 */
public class TraverseActiveRegionsUnitTest extends BaseTest {
    private final static boolean ENFORCE_CONTRACTS = false;
    private final static boolean DEBUG = false;

    @DataProvider(name = "TraversalEngineProvider")
    public Object[][] makeTraversals() {
        final List<Object[]> traversals = new LinkedList<Object[]>();
        traversals.add(new Object[]{new TraverseActiveRegions<>()});
        return traversals.toArray(new Object[][]{});
    }

    private File referenceFile;
    private IndexedFastaSequenceFile reference;
    private SAMSequenceDictionary dictionary;
    private GenomeLocParser genomeLocParser;

    private List<GenomeLoc> intervals;

    private File testBAM;

    @BeforeClass
    private void init() throws IOException {
        //reference = new CachingIndexedFastaSequenceFile(new File("/Users/depristo/Desktop/broadLocal/localData/human_g1k_v37.fasta")); // hg19Reference));
        referenceFile = new File(hg19Reference);
        reference = new CachingIndexedFastaSequenceFile(referenceFile);
        dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);

        // TODO: reads with indels
        // TODO: reads which span many regions
        // TODO: reads which are partially between intervals (in/outside extension)
        // TODO: duplicate reads
        // TODO: read at the end of a contig
        // TODO: reads which are completely outside intervals but within extension
        // TODO: test the extension itself
        // TODO: unmapped reads

        intervals = new ArrayList<GenomeLoc>();
        intervals.add(genomeLocParser.createGenomeLoc("1", 10, 20));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1, 999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        intervals.add(genomeLocParser.createGenomeLoc("2", 1, 100));
        intervals.add(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        intervals = IntervalUtils.sortAndMergeIntervals(genomeLocParser, intervals, IntervalMergingRule.OVERLAPPING_ONLY).toList();

        List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
        reads.add(buildSAMRecord("simple", "1", 100, 200));
        reads.add(buildSAMRecord("overlap_equal", "1", 10, 20));
        reads.add(buildSAMRecord("overlap_unequal", "1", 10, 21));
        reads.add(buildSAMRecord("boundary_equal", "1", 1990, 2009));
        reads.add(buildSAMRecord("boundary_unequal", "1", 1990, 2008));
        reads.add(buildSAMRecord("boundary_1_pre", "1", 1950, 2000));
        reads.add(buildSAMRecord("boundary_1_post", "1", 1999, 2050));
        reads.add(buildSAMRecord("extended_and_np", "1", 990, 1990));
        reads.add(buildSAMRecord("outside_intervals", "1", 5000, 6000));
        reads.add(buildSAMRecord("shard_boundary_1_pre", "1", 16300, 16385));
        reads.add(buildSAMRecord("shard_boundary_1_post", "1", 16384, 16400));
        reads.add(buildSAMRecord("shard_boundary_equal", "1", 16355, 16414));
        reads.add(buildSAMRecord("simple20", "20", 10025, 10075));

        createBAM(reads);
    }

    private void createBAM(List<GATKSAMRecord> reads) throws IOException {
        testBAM = createTempFile("TraverseActiveRegionsUnitTest", ".bam");

        SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reads.get(0).getHeader(), true, testBAM);
        for (GATKSAMRecord read : ReadUtils.sortReadsByCoordinate(reads)) {
            out.addAlignment(read);
        }
        out.close();

        new File(testBAM.getAbsolutePath().replace(".bam", ".bai")).deleteOnExit();
        new File(testBAM.getAbsolutePath() + ".bai").deleteOnExit();
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testAllBasesSeen(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        List<GenomeLoc> activeIntervals = getIsActiveIntervals(t, walker, intervals);
        // Contract: Every genome position in the analysis interval(s) is processed by the walker's isActive() call
        verifyEqualIntervals(intervals, activeIntervals);
    }

    private List<GenomeLoc> getIsActiveIntervals(final TraverseActiveRegions t, DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        List<GenomeLoc> activeIntervals = new ArrayList<GenomeLoc>();
        for (LocusShardDataProvider dataProvider : createDataProviders(t, walker, intervals, testBAM)) {
            t.traverse(walker, dataProvider, 0);
            activeIntervals.addAll(walker.isActiveCalls);
        }

        return activeIntervals;
    }

    @Test (enabled = ENFORCE_CONTRACTS, dataProvider = "TraversalEngineProvider", expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeLow (TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(-0.1);
        getActiveRegions(t, walker, intervals).values();
    }

    @Test (enabled = ENFORCE_CONTRACTS, dataProvider = "TraversalEngineProvider", expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeHigh (TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(1.1);
        getActiveRegions(t, walker, intervals).values();
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testActiveRegionCoverage(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(new GenomeLocSortedSet(genomeLocParser, intervals), true);

        Collection<ActiveRegion> activeRegions = getActiveRegions(t, walker, intervals).values();
        verifyActiveRegionCoverage(intervals, activeRegions);
    }

    private void verifyActiveRegionCoverage(List<GenomeLoc> intervals, Collection<ActiveRegion> activeRegions) {
        List<GenomeLoc> intervalStarts = new ArrayList<GenomeLoc>();
        List<GenomeLoc> intervalStops = new ArrayList<GenomeLoc>();

        for (GenomeLoc interval : intervals) {
            intervalStarts.add(interval.getStartLocation());
            intervalStops.add(interval.getStopLocation());
        }

        Map<GenomeLoc, ActiveRegion> baseRegionMap = new HashMap<GenomeLoc, ActiveRegion>();

        for (ActiveRegion activeRegion : activeRegions) {
            for (GenomeLoc activeLoc : toSingleBaseLocs(activeRegion.getLocation())) {
                // Contract: Regions do not overlap
                Assert.assertFalse(baseRegionMap.containsKey(activeLoc), "Genome location " + activeLoc + " is assigned to more than one region");
                baseRegionMap.put(activeLoc, activeRegion);
            }

            GenomeLoc start = activeRegion.getLocation().getStartLocation();
            if (intervalStarts.contains(start))
                intervalStarts.remove(start);

            GenomeLoc stop = activeRegion.getLocation().getStopLocation();
            if (intervalStops.contains(stop))
                intervalStops.remove(stop);
        }

        for (GenomeLoc baseLoc : toSingleBaseLocs(intervals)) {
            // Contract: Each location in the interval(s) is in exactly one region
            // Contract: The total set of regions exactly matches the analysis interval(s)
            Assert.assertTrue(baseRegionMap.containsKey(baseLoc), "Genome location " + baseLoc + " is not assigned to any region");
            baseRegionMap.remove(baseLoc);
        }

        // Contract: The total set of regions exactly matches the analysis interval(s)
        Assert.assertEquals(baseRegionMap.size(), 0, "Active regions contain base(s) outside of the given intervals");

        // Contract: All explicit interval boundaries must also be region boundaries
        Assert.assertEquals(intervalStarts.size(), 0, "Interval start location does not match an active region start location");
        Assert.assertEquals(intervalStops.size(), 0, "Interval stop location does not match an active region stop location");
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testActiveRegionExtensionOnContig(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        Collection<ActiveRegion> activeRegions = getActiveRegions(t, walker, intervals).values();
        for (ActiveRegion activeRegion : activeRegions) {
            GenomeLoc loc = activeRegion.getExtendedLoc();

            // Contract: active region extensions must stay on the contig
            Assert.assertTrue(loc.getStart() > 0, "Active region extension begins at location " + loc.getStart() + ", past the left end of the contig");
            int refLen = dictionary.getSequence(loc.getContigIndex()).getSequenceLength();
            Assert.assertTrue(loc.getStop() <= refLen, "Active region extension ends at location " + loc.getStop() + ", past the right end of the contig");
        }
    }

    @Test(enabled = true && !DEBUG, dataProvider = "TraversalEngineProvider")
    public void testPrimaryReadMapping(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(new GenomeLocSortedSet(genomeLocParser, intervals),
                EnumSet.of(ActiveRegionReadState.PRIMARY),
                true);

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // extended_and_np: Primary in 1:1-999, Non-Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_equal: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(t, walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal", "extended_and_np");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_unequal", "boundary_1_pre", "boundary_equal", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region);

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testNonPrimaryReadMapping(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(new GenomeLocSortedSet(genomeLocParser, intervals),
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY),
                true);

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // extended_and_np: Primary in 1:1-999, Non-Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_equal: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(t, walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal", "extended_and_np");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testExtendedReadMapping(TraverseActiveRegions t) {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(new GenomeLocSortedSet(genomeLocParser, intervals),
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY, ActiveRegionReadState.EXTENDED),
                true);

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps
        // Contract: Each read has the Extended state in regions where it only overlaps if the region is extended

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_pre: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // boundary_1_post: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // shard_boundary_1_pre: Primary in 1:14908-16384, Non-Primary in 1:16385-16927
        // shard_boundary_1_post: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // shard_boundary_equal: Non-Primary in 1:14908-16384, Primary in 1:16385-16927
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(t, walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));
        verifyReadMapping(region, "simple", "overlap_equal", "overlap_unequal", "extended_and_np");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        verifyReadMapping(region, "boundary_equal", "boundary_unequal", "extended_and_np", "boundary_1_pre", "boundary_1_post");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        verifyReadMapping(region, "shard_boundary_1_pre", "shard_boundary_1_post", "shard_boundary_equal");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        verifyReadMapping(region, "simple20");
    }

    @Test(enabled = true && ! DEBUG, dataProvider = "TraversalEngineProvider")
    public void testUnmappedReads(TraverseActiveRegions t) {
        // TODO
    }

    private void verifyReadMapping(ActiveRegion region, String... reads) {
        Assert.assertNotNull(region, "Region was unexpectedly null");
        final Set<String> regionReads = new HashSet<String>();
        for (SAMRecord read : region.getReads()) {
            Assert.assertFalse(regionReads.contains(read.getReadName()), "Duplicate reads detected in region " + region + " read " + read.getReadName());
            regionReads.add(read.getReadName());
        }

        Collection<String> wantReads = new ArrayList<String>(Arrays.asList(reads));
        for (SAMRecord read : region.getReads()) {
            String regionReadName = read.getReadName();
            Assert.assertTrue(wantReads.contains(regionReadName), "Read " + regionReadName + " incorrectly assigned to active region " + region);
            wantReads.remove(regionReadName);
        }

        Assert.assertTrue(wantReads.isEmpty(), "Reads missing in active region " + region + ", wanted " + (wantReads.isEmpty() ? "" : wantReads.iterator().next()));
    }

    private Map<GenomeLoc, ActiveRegion> getActiveRegions(TraverseActiveRegions t, DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        return getActiveRegions(t, walker, intervals, testBAM);
    }

    private Map<GenomeLoc, ActiveRegion> getActiveRegions(TraverseActiveRegions t, DummyActiveRegionWalker walker, List<GenomeLoc> intervals, final File bam) {
        for (LocusShardDataProvider dataProvider : createDataProviders(t, walker, intervals, bam))
            t.traverse(walker, dataProvider, 0);

        return walker.mappedActiveRegions;
    }

    private Collection<GenomeLoc> toSingleBaseLocs(GenomeLoc interval) {
        List<GenomeLoc> bases = new ArrayList<GenomeLoc>();
        if (interval.size() == 1)
            bases.add(interval);
        else {
            for (int location = interval.getStart(); location <= interval.getStop(); location++)
                bases.add(genomeLocParser.createGenomeLoc(interval.getContig(), location, location));
        }

        return bases;
    }

    private Collection<GenomeLoc> toSingleBaseLocs(List<GenomeLoc> intervals) {
        Set<GenomeLoc> bases = new TreeSet<GenomeLoc>();    // for sorting and uniqueness
        for (GenomeLoc interval : intervals)
            bases.addAll(toSingleBaseLocs(interval));

        return bases;
    }

    private void verifyEqualIntervals(List<GenomeLoc> aIntervals, List<GenomeLoc> bIntervals) {
        Collection<GenomeLoc> aBases = toSingleBaseLocs(aIntervals);
        Collection<GenomeLoc> bBases = toSingleBaseLocs(bIntervals);

        Assert.assertTrue(aBases.size() == bBases.size(), "Interval lists have a differing number of bases: " + aBases.size() + " vs. " + bBases.size());

        Iterator<GenomeLoc> aIter = aBases.iterator();
        Iterator<GenomeLoc> bIter = bBases.iterator();
        while (aIter.hasNext() && bIter.hasNext()) {
            GenomeLoc aLoc = aIter.next();
            GenomeLoc bLoc = bIter.next();
            Assert.assertTrue(aLoc.equals(bLoc), "Interval locations do not match: " + aLoc + " vs. " + bLoc);
        }
    }

    // copied from LocusViewTemplate
    protected GATKSAMRecord buildSAMRecord(String readName, String contig, int alignmentStart, int alignmentEnd) {
        SAMFileHeader header = ArtificialSAMUtils.createDefaultReadGroup(new SAMFileHeader(), "test", "test");
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        GATKSAMRecord record = new GATKSAMRecord(header);

        record.setReadName(readName);
        record.setReferenceIndex(dictionary.getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);

        Cigar cigar = new Cigar();
        int len = alignmentEnd - alignmentStart + 1;
        cigar.add(new CigarElement(len, CigarOperator.M));
        record.setCigar(cigar);
        record.setReadString(new String(new char[len]).replace("\0", "A"));
        record.setBaseQualities(new byte[len]);
        record.setReadGroup(new GATKSAMReadGroupRecord(header.getReadGroup("test")));

        return record;
    }

    private List<LocusShardDataProvider> createDataProviders(TraverseActiveRegions traverseActiveRegions, final Walker walker, List<GenomeLoc> intervals, File bamFile) {
        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);

        Collection<SAMReaderID> samFiles = new ArrayList<SAMReaderID>();
        SAMReaderID readerID = new SAMReaderID(bamFile, new Tags());
        samFiles.add(readerID);

        SAMDataSource dataSource = new SAMDataSource(referenceFile, samFiles, new ThreadAllocation(), null, genomeLocParser,
                false,
                ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                new ArrayList<ReadTransformer>(),
                false, (byte)30, false, true, null, IntervalMergingRule.ALL);

        engine.setReadsDataSource(dataSource);
        final Set<String> samples = ReadUtils.getSAMFileSamples(dataSource.getHeader());

        traverseActiveRegions.initialize(engine, walker);
        List<LocusShardDataProvider> providers = new ArrayList<LocusShardDataProvider>();
        for (Shard shard : dataSource.createShardIteratorOverIntervals(new GenomeLocSortedSet(genomeLocParser, intervals), new ActiveRegionShardBalancer())) {
            for (WindowMaker.WindowMakerIterator window : new WindowMaker(shard, genomeLocParser, dataSource.seek(shard), shard.getGenomeLocs(), samples)) {
                providers.add(new LocusShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>()));
            }
        }

        return providers;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Combinatorial tests to ensure reads are going into the right regions
    //
    // ---------------------------------------------------------------------------------------------------------

    @DataProvider(name = "CombinatorialARTTilingProvider")
    public Object[][] makeCombinatorialARTTilingProvider() {
        final List<Object[]> tests = new LinkedList<Object[]>();

        final List<Integer> starts = Arrays.asList(
                1, // very start of the chromosome
                ArtificialBAMBuilder.BAM_SHARD_SIZE - 100, // right before the shard boundary
                ArtificialBAMBuilder.BAM_SHARD_SIZE + 100 // right after the shard boundary
        );

        final List<EnumSet<ActiveRegionReadState>> allReadStates = Arrays.asList(
                EnumSet.of(ActiveRegionReadState.PRIMARY),
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY),
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY, ActiveRegionReadState.EXTENDED)
        );

        final int maxTests = Integer.MAX_VALUE;
        int nTests = 0;
        for ( final int readLength : Arrays.asList(100) ) {
            for ( final int skips : Arrays.asList(0, 10) ) {
                for ( final int start : starts ) {
                    for ( final int nReadsPerLocus : Arrays.asList(1, 2) ) {
                        for ( final int nLoci : Arrays.asList(1, 1000) ) {
                            final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(reference, nReadsPerLocus, nLoci);
                            bamBuilder.setReadLength(readLength);
                            bamBuilder.setSkipNLoci(skips);
                            bamBuilder.setAlignmentStart(start);
                            for ( EnumSet<ActiveRegionReadState> readStates : allReadStates ) {
                                for ( final GenomeLocSortedSet activeRegions : enumerateActiveRegions(bamBuilder.getAlignmentStart(), bamBuilder.getAlignmentEnd())) {
                                    nTests++;
                                    if ( nTests < maxTests ) // && nTests == 1238 )
                                        tests.add(new Object[]{new TraverseActiveRegions<>(), nTests, activeRegions, readStates, bamBuilder});
                                }
                            }
                        }
                    }
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private Collection<GenomeLocSortedSet> enumerateActiveRegions(final int start, final int stop) {
        // should basically cut up entire region into equal sized chunks, of
        // size 10, 20, 50, 100, etc, alternating skipping pieces so they are inactive
        // Need to make sure we include some edge cases:
        final List<GenomeLocSortedSet> activeRegions = new LinkedList<GenomeLocSortedSet>();

        for ( final int stepSize : Arrays.asList(11, 29, 53, 97) ) {
            for ( final boolean startWithActive : Arrays.asList(true, false) ) {
                activeRegions.add(makeActiveRegionMask(start, stop, stepSize,  startWithActive));
            }
        }

        // active region is the whole interval
        activeRegions.add(new GenomeLocSortedSet(genomeLocParser, genomeLocParser.createGenomeLoc("1", start, stop)));

        // active region extends up to the end of the data, but doesn't include start
        activeRegions.add(new GenomeLocSortedSet(genomeLocParser, genomeLocParser.createGenomeLoc("1", start+10, stop)));

        return activeRegions;
    }

    private GenomeLocSortedSet makeActiveRegionMask(final int start, final int stop, final int stepSize, final boolean startWithActive) {
        final GenomeLocSortedSet active = new GenomeLocSortedSet(genomeLocParser);

        boolean includeRegion = startWithActive;
        for ( int left = start; left < stop; left += stepSize) {
            final int right = left + stepSize;
            final GenomeLoc region = genomeLocParser.createGenomeLoc("1", left, right);
            if ( includeRegion )
                active.add(region);
            includeRegion = ! includeRegion;
        }

        return active;
    }


    @Test(enabled = true && ! DEBUG, dataProvider = "CombinatorialARTTilingProvider")
    public void testARTReadsInActiveRegions(final TraverseActiveRegions<Integer, Integer> traversal, final int id, final GenomeLocSortedSet activeRegions, final EnumSet<ActiveRegionReadState> readStates, final ArtificialBAMBuilder bamBuilder) {
        logger.warn("Running testARTReadsInActiveRegions id=" + id + " locs " + activeRegions + " against bam " + bamBuilder);
        final List<GenomeLoc> intervals = Arrays.asList(
                genomeLocParser.createGenomeLoc("1", bamBuilder.getAlignmentStart(), bamBuilder.getAlignmentEnd())
        );

        final DummyActiveRegionWalker walker = new DummyActiveRegionWalker(activeRegions, false);
        walker.setStates(readStates);

        final Map<GenomeLoc, ActiveRegion> activeRegionsMap = getActiveRegions(traversal, walker, intervals, bamBuilder.makeTemporarilyBAMFile());

        final Set<String> alreadySeenReads = new HashSet<String>(); // for use with the primary / non-primary
        for ( final ActiveRegion region : activeRegionsMap.values() ) {
            final Set<String> readNamesInRegion = readNamesInRegion(region);
            int nReadsExpectedInRegion = 0;
            for ( final GATKSAMRecord read : bamBuilder.makeReads() ) {
                final GenomeLoc readLoc = genomeLocParser.createGenomeLoc(read);

                boolean shouldBeInRegion = readStates.contains(ActiveRegionReadState.EXTENDED)
                        ? region.getExtendedLoc().overlapsP(readLoc)
                        : region.getLocation().overlapsP(readLoc);

                if ( ! readStates.contains(ActiveRegionReadState.NONPRIMARY) ) {
                    if ( alreadySeenReads.contains(read.getReadName()) )
                        shouldBeInRegion = false;
                    else if ( shouldBeInRegion )
                        alreadySeenReads.add(read.getReadName());
                }

                String msg = readNamesInRegion.contains(read.getReadName()) == shouldBeInRegion ? "" : "Region " + region +
                        " failed contains read check: read " + read + " with span " + readLoc + " should be in region is " + shouldBeInRegion + " but I got the opposite";
                Assert.assertEquals(readNamesInRegion.contains(read.getReadName()), shouldBeInRegion, msg);

                nReadsExpectedInRegion += shouldBeInRegion ? 1 : 0;
            }

            Assert.assertEquals(region.size(), nReadsExpectedInRegion, "There are more reads in active region " + region + "than expected");
        }
    }

    private Set<String> readNamesInRegion(final ActiveRegion region) {
        final Set<String> readNames = new LinkedHashSet<String>(region.getReads().size());
        for ( final SAMRecord read : region.getReads() )
            readNames.add(read.getReadName());
        return readNames;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Make sure all insertion reads are properly included in the active regions
    //
    // ---------------------------------------------------------------------------------------------------------

    @Test(dataProvider = "TraversalEngineProvider", enabled = true && ! DEBUG)
    public void ensureAllInsertionReadsAreInActiveRegions(final TraverseActiveRegions<Integer, Integer> traversal) {

        final int readLength = 10;
        final int start = 20;
        final int nReadsPerLocus = 10;
        final int nLoci = 3;

        final ArtificialBAMBuilder bamBuilder = new ArtificialBAMBuilder(reference, nReadsPerLocus, nLoci);
        bamBuilder.setReadLength(readLength);
        bamBuilder.setAlignmentStart(start);

        // note that the position must be +1 as the read's all I cigar puts the end 1 bp before start, leaving it out of the region
        GATKSAMRecord allI = ArtificialSAMUtils.createArtificialRead(bamBuilder.getHeader(),"allI",0,start+1,readLength);
        allI.setCigarString(readLength + "I");
        allI.setReadGroup(new GATKSAMReadGroupRecord(bamBuilder.getHeader().getReadGroups().get(0)));

        bamBuilder.addReads(allI);

        final GenomeLocSortedSet activeRegions = new GenomeLocSortedSet(bamBuilder.getGenomeLocParser());
        activeRegions.add(bamBuilder.getGenomeLocParser().createGenomeLoc("1", 10, 30));
        final List<GenomeLoc> intervals = Arrays.asList(
                genomeLocParser.createGenomeLoc("1", bamBuilder.getAlignmentStart(), bamBuilder.getAlignmentEnd())
        );

        final DummyActiveRegionWalker walker = new DummyActiveRegionWalker(activeRegions, false);

        final Map<GenomeLoc, ActiveRegion> activeRegionsMap = getActiveRegions(traversal, walker, intervals, bamBuilder.makeTemporarilyBAMFile());

        final ActiveRegion region = activeRegionsMap.values().iterator().next();
        int nReadsExpectedInRegion = 0;

        final Set<String> readNamesInRegion = readNamesInRegion(region);
        for ( final GATKSAMRecord read : bamBuilder.makeReads() ) {
            Assert.assertTrue(readNamesInRegion.contains(read.getReadName()),
                    "Region " + region + " should contain read " + read + " with cigar " + read.getCigarString() + " but it wasn't");
            nReadsExpectedInRegion++;
        }

        Assert.assertEquals(region.size(), nReadsExpectedInRegion, "There are more reads in active region " + region + "than expected");
    }
}
