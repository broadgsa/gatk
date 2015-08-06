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

package org.broadinstitute.gatk.utils.refdata;

import htsjdk.samtools.SAMFileHeader;
import org.apache.log4j.Logger;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.codecs.table.TableFeature;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.testng.Assert;
import org.testng.annotations.*;
import java.util.*;
import java.util.List;

public class RefMetaDataTrackerUnitTest {
    final protected static Logger logger = Logger.getLogger(RefMetaDataTrackerUnitTest.class);
    private static SAMFileHeader header;
    private ReferenceContext context;
    private GenomeLocParser genomeLocParser;
    private GenomeLoc locus;
    private final static int START_POS = 10;
    Allele A,C,G,T;
    VariantContext AC_SNP, AG_SNP, AT_SNP;
    TableFeature span10_10, span1_20, span10_20;

    @BeforeClass
    public void beforeClass() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 100);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        locus = genomeLocParser.createGenomeLoc("chr1", START_POS, START_POS);
        context = new ReferenceContext(genomeLocParser, locus, (byte)'A');
        A = Allele.create("A", true);
        C = Allele.create("C");
        G = Allele.create("G");
        T = Allele.create("T");
        AC_SNP = new VariantContextBuilder("x", "chr1", START_POS, START_POS, Arrays.asList(A, C)).make();
        AG_SNP = new VariantContextBuilder("x", "chr1", START_POS, START_POS, Arrays.asList(A, G)).make();
        AT_SNP = new VariantContextBuilder("x", "chr1", START_POS, START_POS, Arrays.asList(A, T)).make();
        span10_10 = makeSpan(10, 10);
        span1_20 = makeSpan(1, 20);
        span10_20 = makeSpan(10, 20);
    }

    @BeforeMethod
    public void reset() {
        RodBinding.resetNameCounter();
    }

    private class MyTest extends BaseTest.TestDataProvider {
        public RODRecordList AValues, BValues;

        private MyTest(Class c, final List<? extends Feature> AValues, final List<? extends Feature> BValues) {
            super(c);
            this.AValues = AValues == null ? null : makeRODRecord("A", AValues);
            this.BValues = BValues == null ? null : makeRODRecord("B", BValues);
        }

        private MyTest(final List<? extends Feature> AValues, final List<? extends Feature> BValues) {
            super(MyTest.class);
            this.AValues = AValues == null ? null : makeRODRecord("A", AValues);
            this.BValues = BValues == null ? null : makeRODRecord("B", BValues);
        }

        @Override
        public String toString() {
            return String.format("A=%s, B=%s", AValues, BValues);
        }

        private final RODRecordList makeRODRecord(String name, List<? extends Feature> features) {
            List<GATKFeature> x = new ArrayList<GATKFeature>();
            for ( Feature f : features )
                x.add(new GATKFeature.TribbleGATKFeature(genomeLocParser, f, name));
            return new RODRecordListImpl(name, x, locus);
        }

        public List<GATKFeature> expected(String name) {
            if ( name.equals("A+B") ) return allValues();
            if ( name.equals("A") ) return expectedAValues();
            if ( name.equals("B") ) return expectedBValues();
            throw new RuntimeException("FAIL");
        }

        public List<GATKFeature> allValues() {
            List<GATKFeature> x = new ArrayList<GATKFeature>();
            x.addAll(expectedAValues());
            x.addAll(expectedBValues());
            return x;
        }

        public List<GATKFeature> expectedAValues() {
            return AValues == null ? Collections.<GATKFeature>emptyList() : AValues;
        }

        public List<GATKFeature> expectedBValues() {
            return BValues == null ? Collections.<GATKFeature>emptyList() : BValues;
        }

        public RefMetaDataTracker makeTracker() {
            List<RODRecordList> x = new ArrayList<RODRecordList>();
            if ( AValues != null ) x.add(AValues);
            if ( BValues != null ) x.add(BValues);
            return new RefMetaDataTracker(x);
        }

        public int nBoundTracks() {
            int n = 0;
            if ( AValues != null ) n++;
            if ( BValues != null ) n++;
            return n;
        }
    }

    private final TableFeature makeSpan(int start, int stop) {
        return new TableFeature(genomeLocParser.createGenomeLoc("chr1", start, stop),
                Collections.<String>emptyList(), Collections.<String>emptyList());
    }

    @DataProvider(name = "tests")
    public Object[][] createTests() {
        new MyTest(null, null);
        new MyTest(Arrays.asList(AC_SNP), null);
        new MyTest(Arrays.asList(AC_SNP, AT_SNP), null);
        new MyTest(Arrays.asList(AC_SNP), Arrays.asList(AG_SNP));
        new MyTest(Arrays.asList(AC_SNP, AT_SNP), Arrays.asList(AG_SNP));
        new MyTest(Arrays.asList(AC_SNP, AT_SNP), Arrays.asList(span10_10));
        new MyTest(Arrays.asList(AC_SNP, AT_SNP), Arrays.asList(span10_10, span10_20));
        new MyTest(Arrays.asList(AC_SNP, AT_SNP), Arrays.asList(span10_10, span10_20, span1_20));

        // for requires starts
        new MyTest(Arrays.asList(span1_20), null);
        new MyTest(Arrays.asList(span10_10, span10_20), null);
        new MyTest(Arrays.asList(span10_10, span10_20, span1_20), null);

        return MyTest.getTests(MyTest.class);
    }

    @Test(enabled = true, dataProvider = "tests")
    public void testRawBindings(MyTest test) {
        logger.warn("Testing " + test + " for number of bound tracks");
        RefMetaDataTracker tracker = test.makeTracker();
        Assert.assertEquals(tracker.getNTracksWithBoundFeatures(), test.nBoundTracks());

        testSimpleBindings("A", tracker, test.AValues);
        testSimpleBindings("B", tracker, test.BValues);
    }

    private <T> void testSimpleBindings(String name, RefMetaDataTracker tracker, RODRecordList expected) {
        List<Feature> asValues = tracker.getValues(Feature.class, name);

        Assert.assertEquals(tracker.hasValues(name), expected != null);
        Assert.assertEquals(asValues.size(), expected == null ? 0 : expected.size());

        if ( expected != null ) {
            for ( GATKFeature e : expected ) {
                boolean foundValue = false;
                for ( Feature f : asValues ) {
                    if ( e.getUnderlyingObject() == f ) foundValue = true;
                }
                Assert.assertTrue(foundValue, "Never found expected value of " + e.getUnderlyingObject() + " bound to " + name + " in " + tracker);
            }
        }
    }

    @Test(enabled = true, dataProvider = "tests")
    public void testGettersAsString(MyTest test) {
        logger.warn("Testing " + test + " for get() methods");
        RefMetaDataTracker tracker = test.makeTracker();

        for ( String name : Arrays.asList("A+B", "A", "B") ) {
            List<Feature> v1 = name.equals("A+B") ? tracker.getValues(Feature.class) : tracker.getValues(Feature.class, name);
            testGetter(name, v1, test.expected(name), true, tracker);

            List<Feature> v2 = name.equals("A+B") ? tracker.getValues(Feature.class, locus) : tracker.getValues(Feature.class, name, locus);
            testGetter(name, v2, startingHere(test.expected(name)), true, tracker);

            Feature v3 = name.equals("A+B") ? tracker.getFirstValue(Feature.class) : tracker.getFirstValue(Feature.class, name);
            testGetter(name, Arrays.asList(v3), test.expected(name), false, tracker);

            Feature v4 = name.equals("A+B") ? tracker.getFirstValue(Feature.class, locus) : tracker.getFirstValue(Feature.class, name, locus);
            testGetter(name, Arrays.asList(v4), startingHere(test.expected(name)), false, tracker);
        }
    }

    @Test(enabled = true, dataProvider = "tests")
    public void testGettersAsRodBindings(MyTest test) {
        logger.warn("Testing " + test + " for get() methods as RodBindings");
        RefMetaDataTracker tracker = test.makeTracker();

        for ( String nameAsString : Arrays.asList("A", "B") ) {
            RodBinding<Feature> binding = new RodBinding<Feature>(Feature.class, nameAsString, "none", "vcf", new Tags());
            List<Feature> v1 = tracker.getValues(binding);
            testGetter(nameAsString, v1, test.expected(nameAsString), true, tracker);

            List<Feature> v2 = tracker.getValues(binding, locus);
            testGetter(nameAsString, v2, startingHere(test.expected(nameAsString)), true, tracker);

            Feature v3 = tracker.getFirstValue(binding);
            testGetter(nameAsString, Arrays.asList(v3), test.expected(nameAsString), false, tracker);

            Feature v4 = tracker.getFirstValue(binding, locus);
            testGetter(nameAsString, Arrays.asList(v4), startingHere(test.expected(nameAsString)), false, tracker);
        }
    }

    @Test(enabled = true, dataProvider = "tests")
    public void testGettersAsListOfRodBindings(MyTest test) {
        logger.warn("Testing " + test + " for get() methods for List<RodBindings>");
        RefMetaDataTracker tracker = test.makeTracker();

        String nameAsString = "A+B";
        RodBinding<Feature> A = new RodBinding<Feature>(Feature.class, "A", "none", "vcf", new Tags());
        RodBinding<Feature> B = new RodBinding<Feature>(Feature.class, "B", "none", "vcf", new Tags());
        List<RodBinding<Feature>> binding = Arrays.asList(A, B);

        List<Feature> v1 = tracker.getValues(binding);
        testGetter(nameAsString, v1, test.expected(nameAsString), true, tracker);

        List<Feature> v2 = tracker.getValues(binding, locus);
        testGetter(nameAsString, v2, startingHere(test.expected(nameAsString)), true, tracker);

        Feature v3 = tracker.getFirstValue(binding);
        testGetter(nameAsString, Arrays.asList(v3), test.expected(nameAsString), false, tracker);

        Feature v4 = tracker.getFirstValue(binding, locus);
        testGetter(nameAsString, Arrays.asList(v4), startingHere(test.expected(nameAsString)), false, tracker);
    }

    private List<GATKFeature> startingHere(List<GATKFeature> l) {
        List<GATKFeature> x = new ArrayList<GATKFeature>();
        for ( GATKFeature f : l ) if ( f.getStart() == locus.getStart() ) x.add(f);
        return x;
    }

    private void testGetter(String name, List<Feature> got, List<GATKFeature> expected, boolean requireExact, RefMetaDataTracker tracker) {
        if ( got.size() == 1 && got.get(0) == null )
            got = Collections.emptyList();

        if ( requireExact )
            Assert.assertEquals(got.size(), expected.size());

        boolean foundAny = false;
        for ( GATKFeature e : expected ) {
            boolean found1 = false;
            for ( Feature got1 : got ) {
                if ( e.getUnderlyingObject() == got1 )
                    found1 = true;
            }
            if ( requireExact )
                Assert.assertTrue(found1, "Never found expected GATKFeature " + e + " bound to " + name + " in " + tracker);
            foundAny = found1 || foundAny;
        }

        if ( ! requireExact && ! expected.isEmpty() )
            Assert.assertTrue(foundAny, "Never found any got values matching one of the expected values bound to " + name + " in " + tracker);
    }
}
