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
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: thibault
 * Date: 11/13/12
 * Time: 2:47 PM
 */
public class TraverseActiveRegionsTest extends BaseTest {

    private class DummyActiveRegionWalker extends ActiveRegionWalker<Integer, Integer> {
        private final double prob;
        public List<GenomeLoc> isActiveCalls = new ArrayList<GenomeLoc>();

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
    private DummyActiveRegionWalker walker;

    @BeforeClass
    private void init() throws FileNotFoundException {
        reference = new CachingIndexedFastaSequenceFile(new File(hg19Reference));
        SAMSequenceDictionary dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);
    }

    @Test
    public void testAllIntervalsSeen() throws Exception {
        List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
        List<GenomeLoc> activeIntervals = new ArrayList<GenomeLoc>();

        GenomeLoc interval = genomeLocParser.createGenomeLoc("1", 1, 1);
        intervals.add(interval);

        LocusShardDataProvider dataProvider = createDataProvider(intervals);

        t.traverse(walker, dataProvider, 0);
        activeIntervals.addAll(walker.isActiveCalls);

        boolean allGenomeLocsSeen = true;
        for (GenomeLoc loc : intervals) {
            boolean thisGenomeLocSeen = false;
            for (GenomeLoc activeLoc : activeIntervals) {
                if (loc.equals(activeLoc)) {
                    thisGenomeLocSeen = true;
                    break;
                }
            }
            if (!thisGenomeLocSeen) {
                allGenomeLocsSeen = false;
                break;
            }
        }

        Assert.assertTrue(allGenomeLocsSeen, "Some intervals missing from activity profile");
    }

    private LocusShardDataProvider createDataProvider(List<GenomeLoc> intervals) {
        walker = new DummyActiveRegionWalker();

        StingSAMIterator iterator = ArtificialSAMUtils.createReadIterator(new ArrayList<SAMRecord>());
        Shard shard = new MockLocusShard(genomeLocParser, intervals);
        WindowMaker windowMaker = new WindowMaker(shard, genomeLocParser,iterator,shard.getGenomeLocs());
        WindowMaker.WindowMakerIterator window = windowMaker.next();

        GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        //engine.setReferenceDataSource(reference);
        engine.setGenomeLocParser(genomeLocParser);
        t.initialize(engine);

        return new LocusShardDataProvider(shard, null, genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>());
    }
}
