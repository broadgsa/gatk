/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.LocusShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.providers.ReadShardDataProvider;
import org.broadinstitute.sting.gatk.datasources.reads.*;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.executive.WindowMaker;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.gatk.traversals.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileState;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.*;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class ReadMetricsUnitTest extends BaseTest {

    @Test
    public void testReadsSeenDoNotOverflowInt() {

        final ReadMetrics metrics = new ReadMetrics();

        final long moreThanMaxInt = ((long)Integer.MAX_VALUE) + 1L;

        for ( long i = 0L; i < moreThanMaxInt; i++ ) {
            metrics.incrementNumReadsSeen();
        }

        Assert.assertEquals(metrics.getNumReadsSeen(), moreThanMaxInt);
        Assert.assertTrue(metrics.getNumReadsSeen() > (long) Integer.MAX_VALUE);

        logger.warn(String.format("%d %d %d", Integer.MAX_VALUE, moreThanMaxInt, Long.MAX_VALUE));
    }


    // Test the accuracy of the read metrics

    private IndexedFastaSequenceFile reference;
    private SAMSequenceDictionary dictionary;
    private SAMFileHeader header;
    private GATKSAMReadGroupRecord readGroup;
    private GenomeLocParser genomeLocParser;
    private File testBAM;

    private static final int numReadsPerContig = 250000;
    private static final List<String> contigs = Arrays.asList("1", "2", "3");

    @BeforeClass
    private void init() throws IOException {
        reference = new CachingIndexedFastaSequenceFile(new File(b37KGReference));
        dictionary = reference.getSequenceDictionary();
        genomeLocParser = new GenomeLocParser(dictionary);
        header = ArtificialSAMUtils.createDefaultReadGroup(new SAMFileHeader(), "test", "test");
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        readGroup = new GATKSAMReadGroupRecord(header.getReadGroup("test"));

        final List<GATKSAMRecord> reads = new ArrayList<>();
        for ( final String contig : contigs ) {
            for ( int i = 1; i <= numReadsPerContig; i++ ) {
                reads.add(buildSAMRecord("read" + contig + "_" + i, contig, i));
            }
        }

        createBAM(reads);
    }

    private void createBAM(final List<GATKSAMRecord> reads) throws IOException {
        testBAM = File.createTempFile("TraverseActiveRegionsUnitTest", ".bam");
        testBAM.deleteOnExit();

        SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reads.get(0).getHeader(), true, testBAM);
        for (GATKSAMRecord read : reads ) {
            out.addAlignment(read);
        }
        out.close();

        new File(testBAM.getAbsolutePath().replace(".bam", ".bai")).deleteOnExit();
        new File(testBAM.getAbsolutePath() + ".bai").deleteOnExit();
    }

    // copied from LocusViewTemplate
    protected GATKSAMRecord buildSAMRecord(final String readName, final String contig, final int alignmentStart) {
        GATKSAMRecord record = new GATKSAMRecord(header);

        record.setReadName(readName);
        record.setReferenceIndex(dictionary.getSequenceIndex(contig));
        record.setAlignmentStart(alignmentStart);

        record.setCigarString("1M");
        record.setReadString("A");
        record.setBaseQualityString("A");
        record.setReadGroup(readGroup);

        return record;
    }

    @Test
    public void testCountsFromReadTraversal() {
        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);

        final Collection<SAMReaderID> samFiles = new ArrayList<>();
        final SAMReaderID readerID = new SAMReaderID(testBAM, new Tags());
        samFiles.add(readerID);

        final SAMDataSource dataSource = new SAMDataSource(samFiles, new ThreadAllocation(), null, genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                new ArrayList<ReadTransformer>(),
                false, (byte)30, false, true);

        engine.setReadsDataSource(dataSource);

        final TraverseReadsNano traverseReadsNano = new TraverseReadsNano(1);
        final DummyReadWalker walker = new DummyReadWalker();
        traverseReadsNano.initialize(engine, walker, null);

        for ( final Shard shard : dataSource.createShardIteratorOverAllReads(new ReadShardBalancer()) ) {
            final ReadShardDataProvider dataProvider = new ReadShardDataProvider(shard, engine.getGenomeLocParser(), dataSource.seek(shard), reference, new ArrayList<ReferenceOrderedDataSource>());
            traverseReadsNano.traverse(walker, dataProvider, 0);
            dataProvider.close();
        }

        Assert.assertEquals(engine.getCumulativeMetrics().getNumReadsSeen(), contigs.size() * numReadsPerContig);
        Assert.assertEquals(engine.getCumulativeMetrics().getNumIterations(), contigs.size() * numReadsPerContig);
    }

    @Test
    public void testCountsFromLocusTraversal() {
        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);

        final Collection<SAMReaderID> samFiles = new ArrayList<>();
        final SAMReaderID readerID = new SAMReaderID(testBAM, new Tags());
        samFiles.add(readerID);

        final SAMDataSource dataSource = new SAMDataSource(samFiles, new ThreadAllocation(), null, genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                new ArrayList<ReadTransformer>(),
                false, (byte)30, false, true);

        engine.setReadsDataSource(dataSource);
        final Set<String> samples = SampleUtils.getSAMFileSamples(dataSource.getHeader());

        final TraverseLociNano traverseLociNano = new TraverseLociNano(1);
        final DummyLocusWalker walker = new DummyLocusWalker();
        traverseLociNano.initialize(engine, walker, null);

        for ( final Shard shard : dataSource.createShardIteratorOverAllReads(new LocusShardBalancer()) ) {
            final WindowMaker windowMaker = new WindowMaker(shard, genomeLocParser, dataSource.seek(shard), shard.getGenomeLocs(), samples);
            for ( WindowMaker.WindowMakerIterator window : windowMaker ) {
                final LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>());
                traverseLociNano.traverse(walker, dataProvider, 0);
                dataProvider.close();
            }
            windowMaker.close();
        }

        //dataSource.close();
        Assert.assertEquals(engine.getCumulativeMetrics().getNumReadsSeen(), contigs.size() * numReadsPerContig);
        Assert.assertEquals(engine.getCumulativeMetrics().getNumIterations(), contigs.size() * numReadsPerContig);
    }

    @Test
    public void testCountsFromActiveRegionTraversal() {
        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);

        final Collection<SAMReaderID> samFiles = new ArrayList<>();
        final SAMReaderID readerID = new SAMReaderID(testBAM, new Tags());
        samFiles.add(readerID);

        final SAMDataSource dataSource = new SAMDataSource(samFiles, new ThreadAllocation(), null, genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                new ArrayList<ReadTransformer>(),
                false, (byte)30, false, true);

        engine.setReadsDataSource(dataSource);
        final Set<String> samples = SampleUtils.getSAMFileSamples(dataSource.getHeader());

        final List<GenomeLoc> intervals = new ArrayList<>(contigs.size());
        for ( final String contig : contigs )
            intervals.add(genomeLocParser.createGenomeLoc(contig, 1, numReadsPerContig));

        final TraverseActiveRegions traverseActiveRegions = new TraverseActiveRegions();
        final DummyActiveRegionWalker walker = new DummyActiveRegionWalker();
        traverseActiveRegions.initialize(engine, walker, null);

        for ( final Shard shard : dataSource.createShardIteratorOverIntervals(new GenomeLocSortedSet(genomeLocParser, intervals), new ActiveRegionShardBalancer()) ) {
            final WindowMaker windowMaker = new WindowMaker(shard, genomeLocParser, dataSource.seek(shard), shard.getGenomeLocs(), samples);
            for ( WindowMaker.WindowMakerIterator window : windowMaker ) {
                final LocusShardDataProvider dataProvider = new LocusShardDataProvider(shard, shard.getReadProperties(), genomeLocParser, window.getLocus(), window, reference, new ArrayList<ReferenceOrderedDataSource>());
                traverseActiveRegions.traverse(walker, dataProvider, 0);
                dataProvider.close();
            }
            windowMaker.close();
        }

        Assert.assertEquals(engine.getCumulativeMetrics().getNumReadsSeen(), contigs.size() * numReadsPerContig);
        Assert.assertEquals(engine.getCumulativeMetrics().getNumIterations(), contigs.size() * numReadsPerContig);
    }

    @Test
    public void testFilteredCounts() {
        final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
        engine.setGenomeLocParser(genomeLocParser);

        final Collection<SAMReaderID> samFiles = new ArrayList<>();
        final SAMReaderID readerID = new SAMReaderID(testBAM, new Tags());
        samFiles.add(readerID);

        final List<ReadFilter> filters = new ArrayList<>();
        filters.add(new EveryTenthReadFilter());

        final SAMDataSource dataSource = new SAMDataSource(samFiles, new ThreadAllocation(), null, genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.STRICT,
                null,
                null,
                new ValidationExclusion(),
                filters,
                new ArrayList<ReadTransformer>(),
                false, (byte)30, false, true);

        engine.setReadsDataSource(dataSource);

        final TraverseReadsNano traverseReadsNano = new TraverseReadsNano(1);
        final DummyReadWalker walker = new DummyReadWalker();
        traverseReadsNano.initialize(engine, walker, null);

        for ( final Shard shard : dataSource.createShardIteratorOverAllReads(new ReadShardBalancer()) ) {
            final ReadShardDataProvider dataProvider = new ReadShardDataProvider(shard, engine.getGenomeLocParser(), dataSource.seek(shard), reference, new ArrayList<ReferenceOrderedDataSource>());
            traverseReadsNano.traverse(walker, dataProvider, 0);
            dataProvider.close();
        }

        Assert.assertEquals((long)engine.getCumulativeMetrics().getCountsByFilter().get(EveryTenthReadFilter.class.getSimpleName()), contigs.size() * numReadsPerContig / 10);
    }

    class DummyLocusWalker extends LocusWalker<Integer, Integer> {
        @Override
        public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
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

    class DummyReadWalker extends ReadWalker<Integer, Integer> {
        @Override
        public Integer map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
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

    class DummyActiveRegionWalker extends ActiveRegionWalker<Integer, Integer> {
        @Override
        public ActivityProfileState isActive(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
            return new ActivityProfileState(ref.getLocus(), 0.0);
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

    private final class EveryTenthReadFilter extends ReadFilter {

        private int myCounter = 0;

        @Override
        public boolean filterOut(final SAMRecord record) {
            if ( ++myCounter == 10 ) {
                myCounter = 0;
                return true;
            }

            return false;
        }
    }
}