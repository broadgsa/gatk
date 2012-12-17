package org.broadinstitute.sting.gatk.traversals;

import com.google.java.contract.PreconditionError;
import net.sf.samtools.*;
import org.broadinstitute.sting.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.reads.MockLocusShard;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.executive.WindowMaker;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;


import java.io.File;
import java.io.FileNotFoundException;
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
public class TraverseActiveRegionsTest extends BaseTest {

    private class DummyActiveRegionWalker extends ActiveRegionWalker<Integer, Integer> {
        private final double prob;
        private EnumSet<ActiveRegionReadState> states = super.desiredReadStates();

        protected List<GenomeLoc> isActiveCalls = new ArrayList<GenomeLoc>();
        protected Map<GenomeLoc, ActiveRegion> mappedActiveRegions = new HashMap<GenomeLoc, ActiveRegion>();

        public DummyActiveRegionWalker() {
            this.prob = 1.0;
        }

        public DummyActiveRegionWalker(double constProb) {
            this.prob = constProb;
        }

        public DummyActiveRegionWalker(EnumSet<ActiveRegionReadState> wantStates) {
            this.prob = 1.0;
            this.states = wantStates;
        }

        @Override
        public EnumSet<ActiveRegionReadState> desiredReadStates() {
            return states;
        }

        @Override
        public ActivityProfileResult isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            isActiveCalls.add(ref.getLocus());
            return new ActivityProfileResult(ref.getLocus(), prob);
        }

        @Override
        public Integer map(ActiveRegion activeRegion, RefMetaDataTracker metaDataTracker) {
            mappedActiveRegions.put(activeRegion.getLocation(), activeRegion);
            return 0;
        }

        @Override
        public Integer reduceInit() {
            return 0;
        }

        @Override
        public Integer reduce(Integer value, Integer sum) {
            return 0;
        }
    }

    private final TraverseActiveRegions<Integer, Integer> t = new TraverseActiveRegions<Integer, Integer>();

    private IndexedFastaSequenceFile reference;
    private SAMSequenceDictionary dictionary;
    private GenomeLocParser genomeLocParser;

    private List<GenomeLoc> intervals;
    private List<GATKSAMRecord> reads;

    @BeforeClass
    private void init() throws FileNotFoundException {
        reference = new CachingIndexedFastaSequenceFile(new File(hg19Reference));
        dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);

        // TODO: test shard boundaries

        intervals = new ArrayList<GenomeLoc>();
        intervals.add(genomeLocParser.createGenomeLoc("1", 10, 20));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1, 999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 2000, 2999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 10000, 20000));
        intervals.add(genomeLocParser.createGenomeLoc("1", 249250600, 249250621));
        intervals.add(genomeLocParser.createGenomeLoc("2", 1, 100));
        intervals.add(genomeLocParser.createGenomeLoc("20", 10000, 10100));
        intervals = IntervalUtils.sortAndMergeIntervals(genomeLocParser, intervals, IntervalMergingRule.OVERLAPPING_ONLY).toList();

        reads = new ArrayList<GATKSAMRecord>();
        reads.add(buildSAMRecord("simple", "1", 100, 200));
        reads.add(buildSAMRecord("overlap_equal", "1", 10, 20));
        reads.add(buildSAMRecord("overlap_unequal", "1", 10, 21));
        reads.add(buildSAMRecord("boundary_equal", "1", 1990, 2009));
        reads.add(buildSAMRecord("boundary_unequal", "1", 1990, 2008));
        reads.add(buildSAMRecord("extended_and_np", "1", 990, 1990));
        reads.add(buildSAMRecord("outside_intervals", "1", 5000, 6000));
        reads.add(buildSAMRecord("end_of_chr1", "1", 249250600, 249250700));
        reads.add(buildSAMRecord("simple20", "20", 10025, 10075));

        // required by LocusIteratorByState, and I prefer to list them in test case order above
        ReadUtils.sortReadsByCoordinate(reads);
    }

    @Test
    public void testAllBasesSeen() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        List<GenomeLoc> activeIntervals = getIsActiveIntervals(walker, intervals);
        // Contract: Every genome position in the analysis interval(s) is processed by the walker's isActive() call
        verifyEqualIntervals(intervals, activeIntervals);
    }

    private List<GenomeLoc> getIsActiveIntervals(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        List<GenomeLoc> activeIntervals = new ArrayList<GenomeLoc>();
        for (LocusShardDataProvider dataProvider : createDataProviders(intervals)) {
            t.traverse(walker, dataProvider, 0);
            activeIntervals.addAll(walker.isActiveCalls);
        }

        return activeIntervals;
    }

    @Test (expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeLow () {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(-0.1);
        getActiveRegions(walker, intervals).values();
    }

    @Test (expectedExceptions = PreconditionError.class)
    public void testIsActiveRangeHigh () {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(1.1);
        getActiveRegions(walker, intervals).values();
    }

    @Test
    public void testActiveRegionCoverage() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        Collection<ActiveRegion> activeRegions = getActiveRegions(walker, intervals).values();
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

    @Test
    public void testPrimaryReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // end_of_chr1: Primary in 1:249250600-249250621
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));

        getRead(region, "simple");
        getRead(region, "overlap_equal");
        getRead(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        getRead(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        getRead(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 249250600, 249250621));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        getRead(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        getRead(region, "simple20");
    }

    @Test
    public void testNonPrimaryReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY));

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // end_of_chr1: Primary in 1:249250600-249250621
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));

        getRead(region, "simple");
        getRead(region, "overlap_equal");
        getRead(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        getRead(region, "boundary_equal");
        getRead(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        getRead(region, "boundary_equal");
        getRead(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 249250600, 249250621));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        getRead(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        getRead(region, "simple20");
    }

    @Test
    public void testExtendedReadMapping() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker(
                EnumSet.of(ActiveRegionReadState.PRIMARY, ActiveRegionReadState.NONPRIMARY, ActiveRegionReadState.EXTENDED));

        // Contract: Each read has the Primary state in a single region (or none)
        // This is the region of maximum overlap for the read (earlier if tied)

        // Contract: Each read has the Non-Primary state in all other regions it overlaps
        // Contract: Each read has the Extended state in regions where it only overlaps if the region is extended

        // simple: Primary in 1:1-999
        // overlap_equal: Primary in 1:1-999
        // overlap_unequal: Primary in 1:1-999
        // boundary_equal: Non-Primary in 1:1000-1999, Primary in 1:2000-2999
        // boundary_unequal: Primary in 1:1000-1999, Non-Primary in 1:2000-2999
        // extended_and_np: Non-Primary in 1:1-999, Primary in 1:1000-1999, Extended in 1:2000-2999
        // outside_intervals: none
        // end_of_chr1: Primary in 1:249250600-249250621
        // simple20: Primary in 20:10000-10100

        Map<GenomeLoc, ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        ActiveRegion region;

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1, 999));

        getRead(region, "simple");
        getRead(region, "overlap_equal");
        getRead(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 1000, 1999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        getRead(region, "boundary_equal");
        getRead(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 2000, 2999));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        getRead(region, "boundary_equal");
        getRead(region, "boundary_unequal");
        getRead(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("1", 249250600, 249250621));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        getRead(region, "end_of_chr1");
        verifyReadNotPlaced(region, "simple20");

        region = activeRegions.get(genomeLocParser.createGenomeLoc("20", 10000, 10100));

        verifyReadNotPlaced(region, "simple");
        verifyReadNotPlaced(region, "overlap_equal");
        verifyReadNotPlaced(region, "overlap_unequal");
        verifyReadNotPlaced(region, "boundary_equal");
        verifyReadNotPlaced(region, "boundary_unequal");
        verifyReadNotPlaced(region, "extended_and_np");
        verifyReadNotPlaced(region, "outside_intervals");
        verifyReadNotPlaced(region, "end_of_chr1");
        getRead(region, "simple20");
    }

    @Test
    public void testUnmappedReads() {
        // TODO
    }

    private void verifyReadNotPlaced(ActiveRegion region, String readName) {
        for (SAMRecord read : region.getReads()) {
            if (read.getReadName().equals(readName))
                Assert.fail("Read " + readName + " found in active region " + region);
        }
    }

    private SAMRecord getRead(ActiveRegion region, String readName) {
        for (SAMRecord read : region.getReads()) {
            if (read.getReadName().equals(readName))
                return read;
        }

        Assert.fail("Read " + readName + " not assigned to active region " + region);
        return null;
    }

    private Map<GenomeLoc, ActiveRegion> getActiveRegions(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        for (LocusShardDataProvider dataProvider : createDataProviders(intervals))
            t.traverse(walker, dataProvider, 0);

        t.endTraversal(walker, 0);

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
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        GATKSAMRecord record = new GATKSAMRecord(header);

        record.setReadName(readName);
        record.setReferenceIndex(dictionary.getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);

        Cigar cigar = new Cigar();
        int len = alignmentEnd - alignmentStart + 1;
        cigar.add(new CigarElement(len, CigarOperator.M));
        record.setCigar(cigar);
        record.setReadBases(new byte[len]);
        record.setBaseQualities(new byte[len]);

        return record;
    }

    private List<LocusShardDataProvider> createDataProviders(List<GenomeLoc> intervals) {
        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);
        t.initialize(engine);

        StingSAMIterator iterator = ArtificialSAMUtils.createReadIterator(new ArrayList<SAMRecord>(reads));
        Shard shard = new MockLocusShard(genomeLocParser, intervals);

        List<LocusShardDataProvider> providers = new ArrayList<LocusShardDataProvider>();
        for (WindowMaker.WindowMakerIterator window : new WindowMaker(shard, genomeLocParser, iterator, shard.getGenomeLocs())) {
            providers.add(new LocusShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>()));
        }

        return providers;
    }
}
