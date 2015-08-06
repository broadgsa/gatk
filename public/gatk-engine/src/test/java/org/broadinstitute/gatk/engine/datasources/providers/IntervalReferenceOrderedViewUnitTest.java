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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.SimpleFeature;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.refdata.RODRecordListImpl;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * @author depristo
 */
public class IntervalReferenceOrderedViewUnitTest extends BaseTest {
    private static int startingChr = 1;
    private static int endingChr = 2;
    private static int readCount = 100;
    private static int DEFAULT_READ_LENGTH = ArtificialSAMUtils.DEFAULT_READ_LENGTH;
    private static String contig;
    private static SAMFileHeader header;

    private GenomeLocParser genomeLocParser;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader((endingChr - startingChr) + 1, startingChr, readCount + DEFAULT_READ_LENGTH);
        contig = header.getSequence(0).getSequenceName();
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());

        initializeTests();
    }

    private class CompareFeatures implements Comparator<Feature> {
        @Override
        public int compare(Feature o1, Feature o2) {
            return genomeLocParser.createGenomeLoc(o1).compareTo(genomeLocParser.createGenomeLoc(o2));
        }
    }

    private class ReadMetaDataTrackerRODStreamTest extends TestDataProvider {
        final List<Feature> allFeatures;
        final List<GenomeLoc> intervals;

        public ReadMetaDataTrackerRODStreamTest(final List<Feature> allFeatures, final GenomeLoc interval) {
            this(allFeatures, Collections.singletonList(interval));
        }

        public ReadMetaDataTrackerRODStreamTest(final List<Feature> allFeatures, final List<GenomeLoc> intervals) {
            super(ReadMetaDataTrackerRODStreamTest.class);
            this.allFeatures = new ArrayList<Feature>(allFeatures);
            Collections.sort(this.allFeatures, new CompareFeatures());
            this.intervals = new ArrayList<GenomeLoc>(intervals);
            Collections.sort(this.intervals);
            setName(String.format("%s nFeatures %d intervals %s", getClass().getSimpleName(), allFeatures.size(),
                    intervals.size() == 1 ? intervals.get(0) : "size " + intervals.size()));
        }

        public PeekableIterator<RODRecordList> getIterator(final String name) {
            return new PeekableIterator<RODRecordList>(new TribbleIteratorFromCollection(name, genomeLocParser, allFeatures));
        }

        public Set<Feature> getExpectedOverlaps(final GenomeLoc interval) {
            final Set<Feature> overlapping = new HashSet<Feature>();
            for ( final Feature f : allFeatures )
                if ( genomeLocParser.createGenomeLoc(f).overlapsP(interval) )
                    overlapping.add(f);
            return overlapping;
        }
    }

    public void initializeTests() {
        final List<Feature> handPickedFeatures = new ArrayList<Feature>();

        handPickedFeatures.add(new SimpleFeature(contig, 1, 1));
        handPickedFeatures.add(new SimpleFeature(contig, 2, 5));
        handPickedFeatures.add(new SimpleFeature(contig, 4, 4));
        handPickedFeatures.add(new SimpleFeature(contig, 6, 6));
        handPickedFeatures.add(new SimpleFeature(contig, 9, 10));
        handPickedFeatures.add(new SimpleFeature(contig, 10, 10));
        handPickedFeatures.add(new SimpleFeature(contig, 10, 11));
        handPickedFeatures.add(new SimpleFeature(contig, 13, 20));

        createTestsForFeatures(handPickedFeatures);

        // test in the present of a large spanning element
        {
            List<Feature> oneLargeSpan = new ArrayList<Feature>(handPickedFeatures);
            oneLargeSpan.add(new SimpleFeature(contig, 1, 30));
            createTestsForFeatures(oneLargeSpan);
        }

        // test in the presence of a partially spanning element
        {
            List<Feature> partialSpanStart = new ArrayList<Feature>(handPickedFeatures);
            partialSpanStart.add(new SimpleFeature(contig, 1, 6));
            createTestsForFeatures(partialSpanStart);
        }

        // test in the presence of a partially spanning element at the end
        {
            List<Feature> partialSpanEnd = new ArrayList<Feature>(handPickedFeatures);
            partialSpanEnd.add(new SimpleFeature(contig, 10, 30));
            createTestsForFeatures(partialSpanEnd);
        }

        // no data at all
        final GenomeLoc loc = genomeLocParser.createGenomeLoc(contig, 5, 5);
        new ReadMetaDataTrackerRODStreamTest(Collections.<Feature>emptyList(), loc);
    }

    // --------------------------------------------------------------------------------
    //
    // tests for the lower level IntervalOverlappingRODsFromStream
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "ReadMetaDataTrackerRODStreamTest")
    public Object[][] createReadMetaDataTrackerRODStreamTest() {
        return ReadMetaDataTrackerRODStreamTest.getTests(ReadMetaDataTrackerRODStreamTest.class);
    }

    private GenomeLoc span(final List<GenomeLoc> features) {
        int featuresStart = 1; for ( final GenomeLoc f : features ) featuresStart = Math.min(featuresStart, f.getStart());
        int featuresStop = 1; for ( final GenomeLoc f : features ) featuresStop = Math.max(featuresStop, f.getStop());
        return genomeLocParser.createGenomeLoc(contig, featuresStart, featuresStop);
    }

    private void createTestsForFeatures(final List<Feature> features) {
        int featuresStart = 1; for ( final Feature f : features ) featuresStart = Math.min(featuresStart, f.getStart());
        int featuresStop = 1; for ( final Feature f : features ) featuresStop = Math.max(featuresStop, f.getEnd());

        for ( final int size : Arrays.asList(1, 5, 10, 100) ) {
            final List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>();
            // regularly spaced
            for ( int start = featuresStart; start < featuresStop; start++) {
                final GenomeLoc loc = genomeLocParser.createGenomeLoc(contig, start, start + size - 1);
                allIntervals.add(loc);
                new ReadMetaDataTrackerRODStreamTest(features, loc);
            }

            // starting and stopping at every feature
            for ( final Feature f : features ) {
                // just at the feature
                allIntervals.add(genomeLocParser.createGenomeLoc(contig, f.getStart(), f.getEnd()));
                new ReadMetaDataTrackerRODStreamTest(features, allIntervals.get(allIntervals.size() - 1));

                // up to end
                allIntervals.add(genomeLocParser.createGenomeLoc(contig, f.getStart() - 1, f.getEnd()));
                new ReadMetaDataTrackerRODStreamTest(features, allIntervals.get(allIntervals.size() - 1));

                // missing by 1
                allIntervals.add(genomeLocParser.createGenomeLoc(contig, f.getStart() + 1, f.getEnd() + 1));
                new ReadMetaDataTrackerRODStreamTest(features, allIntervals.get(allIntervals.size() - 1));

                // just spanning
                allIntervals.add(genomeLocParser.createGenomeLoc(contig, f.getStart() - 1, f.getEnd() + 1));
                new ReadMetaDataTrackerRODStreamTest(features, allIntervals.get(allIntervals.size() - 1));
            }

            new ReadMetaDataTrackerRODStreamTest(features, allIntervals);
        }
    }

    @Test(enabled = true, dataProvider = "ReadMetaDataTrackerRODStreamTest")
    public void runReadMetaDataTrackerRODStreamTest_singleQuery(final ReadMetaDataTrackerRODStreamTest data) {
        if ( data.intervals.size() == 1 ) {
            final String name = "testName";
            final PeekableIterator<RODRecordList> iterator = data.getIterator(name);
            final IntervalOverlappingRODsFromStream stream = new IntervalOverlappingRODsFromStream(name, iterator);
            testRODStream(data, stream, Collections.singletonList(data.intervals.get(0)));
        }
    }

    @Test(enabled = true, dataProvider = "ReadMetaDataTrackerRODStreamTest", dependsOnMethods = "runReadMetaDataTrackerRODStreamTest_singleQuery")
    public void runReadMetaDataTrackerRODStreamTest_multipleQueries(final ReadMetaDataTrackerRODStreamTest data) {
        if ( data.intervals.size() > 1 ) {
            final String name = "testName";
            final PeekableIterator<RODRecordList> iterator = data.getIterator(name);
            final IntervalOverlappingRODsFromStream stream = new IntervalOverlappingRODsFromStream(name, iterator);
            testRODStream(data, stream, data.intervals);
        }
    }

    private void testRODStream(final ReadMetaDataTrackerRODStreamTest test, final IntervalOverlappingRODsFromStream stream, final List<GenomeLoc> intervals) {
        for ( final GenomeLoc interval : intervals ) {
            final RODRecordList query = stream.getOverlapping(interval);
            final HashSet<Feature> queryFeatures = new HashSet<Feature>();
            for ( final GATKFeature f : query ) queryFeatures.add((Feature)f.getUnderlyingObject());
            final Set<Feature> overlaps = test.getExpectedOverlaps(interval);

            Assert.assertEquals(queryFeatures.size(), overlaps.size(), "IntervalOverlappingRODsFromStream didn't return the expected set of overlapping features." +
                    " Expected size = " + overlaps.size() + " but saw " + queryFeatures.size());

            BaseTest.assertEqualsSet(queryFeatures, overlaps, "IntervalOverlappingRODsFromStream didn't return the expected set of overlapping features." +
                    " Expected = " + Utils.join(",", overlaps) + " but saw " + Utils.join(",", queryFeatures));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // tests for the higher level tracker itself
    //
    // --------------------------------------------------------------------------------

    @DataProvider(name = "ReadMetaDataTrackerTests")
    public Object[][] createTrackerTests() {
        List<Object[]> tests = new ArrayList<Object[]>();

        final Object[][] singleTests = ReadMetaDataTrackerRODStreamTest.getTests(ReadMetaDataTrackerRODStreamTest.class);
        final List<ReadMetaDataTrackerRODStreamTest> multiSiteTests = new ArrayList<ReadMetaDataTrackerRODStreamTest>();
        for ( final Object[] singleTest : singleTests ) {
            if ( ((ReadMetaDataTrackerRODStreamTest)singleTest[0]).intervals.size() > 1 )
                multiSiteTests.add((ReadMetaDataTrackerRODStreamTest)singleTest[0]);
        }

        for ( final boolean testStateless : Arrays.asList(true, false) ) {
            // all pairwise tests
            for ( List<ReadMetaDataTrackerRODStreamTest> singleTest : Utils.makePermutations(multiSiteTests, 2, false)) {
                tests.add(new Object[]{singleTest, testStateless});
            }

            // all 3 way pairwise tests
            //for ( List<ReadMetaDataTrackerRODStreamTest> singleTest : Utils.makePermutations(multiSiteTests, 3, false)) {
            //    tests.add(new Object[]{singleTest, testStateless});
            //}
        }

        logger.warn("Creating " + tests.size() + " tests for ReadMetaDataTrackerTests");
        return tests.toArray(new Object[][]{});
    }

    @Test(enabled = true, dataProvider = "ReadMetaDataTrackerTests", dependsOnMethods = "runReadMetaDataTrackerRODStreamTest_multipleQueries")
    public void runReadMetaDataTrackerTest(final List<ReadMetaDataTrackerRODStreamTest> RODs, final boolean testStateless) {
        final List<String> names = new ArrayList<String>();
        final List<PeekableIterator<RODRecordList>> iterators = new ArrayList<PeekableIterator<RODRecordList>>();
        final List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
        final List<RodBinding<Feature>> rodBindings = new ArrayList<RodBinding<Feature>>();

        for ( int i = 0; i < RODs.size(); i++ ) {
            final RodBinding<Feature> rodBinding = new RodBinding<Feature>(Feature.class, "name"+i);
            rodBindings.add(rodBinding);
            final String name = rodBinding.getName();
            names.add(name);
            iterators.add(RODs.get(i).getIterator(name));
            intervals.addAll(RODs.get(i).intervals);
        }

        Collections.sort(intervals);
        final GenomeLoc span = span(intervals);
        final IntervalReferenceOrderedView view = new IntervalReferenceOrderedView(genomeLocParser, span, names, iterators);

        if ( testStateless ) {
            // test each tracker is well formed, as each is created
            for ( final GenomeLoc interval : intervals ) {
                final RefMetaDataTracker tracker = view.getReferenceOrderedDataForInterval(interval);
                testMetaDataTrackerBindings(tracker, interval, RODs, rodBindings);
            }
        } else {
            // tests all trackers are correct after reading them into an array
            // this checks that the trackers are be safely stored away and analyzed later (critical for nano-scheduling)
            final List<RefMetaDataTracker> trackers = new ArrayList<RefMetaDataTracker>();
            for ( final GenomeLoc interval : intervals ) {
                final RefMetaDataTracker tracker = view.getReferenceOrderedDataForInterval(interval);
                trackers.add(tracker);
            }

            for ( int i = 0; i < trackers.size(); i++) {
                testMetaDataTrackerBindings(trackers.get(i), intervals.get(i), RODs, rodBindings);
            }
        }
    }

    private void testMetaDataTrackerBindings(final RefMetaDataTracker tracker,
                                             final GenomeLoc interval,
                                             final List<ReadMetaDataTrackerRODStreamTest> RODs,
                                             final List<RodBinding<Feature>> rodBindings) {
        for ( int i = 0; i < RODs.size(); i++ ) {
            final ReadMetaDataTrackerRODStreamTest test = RODs.get(i);
            final List<Feature> queryFeaturesList = tracker.getValues(rodBindings.get(i));
            final Set<Feature> queryFeatures = new HashSet<Feature>(queryFeaturesList);
            final Set<Feature> overlaps = test.getExpectedOverlaps(interval);

            Assert.assertEquals(queryFeatures.size(), overlaps.size(), "IntervalOverlappingRODsFromStream didn't return the expected set of overlapping features." +
                    " Expected size = " + overlaps.size() + " but saw " + queryFeatures.size());

            BaseTest.assertEqualsSet(queryFeatures, overlaps, "IntervalOverlappingRODsFromStream didn't return the expected set of overlapping features." +
                    " Expected = " + Utils.join(",", overlaps) + " but saw " + Utils.join(",", queryFeatures));
        }
    }

    static class TribbleIteratorFromCollection implements Iterator<RODRecordList> {
        // current location
        private final String name;
        final Queue<GATKFeature> gatkFeatures;

        public TribbleIteratorFromCollection(final String name, final GenomeLocParser genomeLocParser, final List<Feature> features) {
            this.name = name;

            this.gatkFeatures = new LinkedList<GATKFeature>();
            for ( final Feature f : features )
                gatkFeatures.add(new GATKFeature.TribbleGATKFeature(genomeLocParser, f, name));
        }

        @Override
        public boolean hasNext() {
            return ! gatkFeatures.isEmpty();
        }

        @Override
        public RODRecordList next() {
            final GATKFeature first = gatkFeatures.poll();
            final Collection<GATKFeature> myFeatures = new LinkedList<GATKFeature>();
            myFeatures.add(first);
            while ( gatkFeatures.peek() != null && gatkFeatures.peek().getLocation().getStart() == first.getStart() )
                myFeatures.add(gatkFeatures.poll());

            GenomeLoc loc = first.getLocation();
            for ( final GATKFeature feature : myFeatures )
                loc = loc.merge(feature.getLocation());

            return new RODRecordListImpl(name, myFeatures, loc); // is this safe?
        }

        @Override public void remove() { throw new IllegalStateException("GRRR"); }
    }
}


