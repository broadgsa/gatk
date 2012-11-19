package org.broadinstitute.sting.gatk.traversals;

import org.testng.Assert;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
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
        public List<GenomeLoc> isActiveCalls = new ArrayList<GenomeLoc>();
        public List<ActiveRegion> mappedActiveRegions = new ArrayList<ActiveRegion>();

        public DummyActiveRegionWalker() {
            this.prob = 1.0;
        }

        @Override
        public ActivityProfileResult isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            isActiveCalls.add(ref.getLocus());
            return new ActivityProfileResult(ref.getLocus(), prob);
        }

        @Override
        public Integer map(ActiveRegion activeRegion, RefMetaDataTracker metaDataTracker) {
            mappedActiveRegions.add(activeRegion);
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
    private GenomeLocParser genomeLocParser;

    @BeforeClass
    private void init() throws FileNotFoundException {
        reference = new CachingIndexedFastaSequenceFile(new File(hg19Reference));
        SAMSequenceDictionary dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);
    }

    @Test
    public void testAllBasesSeen() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();
        List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();

        intervals.add(genomeLocParser.createGenomeLoc("1", 1, 1));
        List<GenomeLoc> activeIntervals = getIsActiveIntervals(walker, intervals);
        // Contract: Every genome position in the analysis interval(s) is processed by the walker's isActive() call
        verifyEqualIntervals(intervals, activeIntervals);

        intervals.add(genomeLocParser.createGenomeLoc("1", 10, 20));
        activeIntervals = getIsActiveIntervals(walker, intervals);
        // Contract: Every genome position in the analysis interval(s) is processed by the walker's isActive() call
        verifyEqualIntervals(intervals, activeIntervals);

        // TODO: more tests and edge cases
    }

    private List<GenomeLoc> getIsActiveIntervals(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        List<GenomeLoc> activeIntervals = new ArrayList<GenomeLoc>();
        for (LocusShardDataProvider dataProvider : createDataProviders(intervals)) {
            t.traverse(walker, dataProvider, 0);
            activeIntervals.addAll(walker.isActiveCalls);
        }

        return activeIntervals;
    }

    @Test
    public void testActiveRegionCoverage() {
        DummyActiveRegionWalker walker = new DummyActiveRegionWalker();
        List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();

        intervals.add(genomeLocParser.createGenomeLoc("1", 1, 999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 1000, 1999));
        intervals.add(genomeLocParser.createGenomeLoc("1", 2000, 2999));

        List<ActiveRegion> activeRegions = getActiveRegions(walker, intervals);
        verifyActiveRegionCoverage(intervals, activeRegions);

        // TODO: more tests and edge cases
    }

    private void verifyActiveRegionCoverage(List<GenomeLoc> intervals, List<ActiveRegion> activeRegions) {
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

    private List<ActiveRegion> getActiveRegions(DummyActiveRegionWalker walker, List<GenomeLoc> intervals) {
        for (LocusShardDataProvider dataProvider : createDataProviders(intervals))
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

    private List<LocusShardDataProvider> createDataProviders(List<GenomeLoc> intervals) {
        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);
        t.initialize(engine);

        StingSAMIterator iterator = ArtificialSAMUtils.createReadIterator(new ArrayList<SAMRecord>());
        Shard shard = new MockLocusShard(genomeLocParser, intervals);

        List<LocusShardDataProvider> providers = new ArrayList<LocusShardDataProvider>();
        for (WindowMaker.WindowMakerIterator window : new WindowMaker(shard, genomeLocParser, iterator, shard.getGenomeLocs())) {
            providers.add(new LocusShardDataProvider(shard, null, genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>()));
        }

        return providers;
    }
}
