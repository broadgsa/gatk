/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package org.broadinstitute.sting.utils.variantcontext;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class VariantContextUtilsUnitTest extends BaseTest {
    Allele Aref, T, C, delRef, Cref, ATC, ATCATC;
    private GenomeLocParser genomeLocParser;

    @BeforeSuite
    public void setup() {
        final File referenceFile = new File(b37KGReference);
        try {
            IndexedFastaSequenceFile seq = new CachingIndexedFastaSequenceFile(referenceFile);
            genomeLocParser = new GenomeLocParser(seq);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }

        // alleles
        Aref = Allele.create("A", true);
        Cref = Allele.create("C", true);
        delRef = Allele.create("-", true);
        T = Allele.create("T");
        C = Allele.create("C");
        ATC = Allele.create("ATC");
        ATCATC = Allele.create("ATCATC");
    }

    private Genotype makeG(String sample, Allele a1, Allele a2) {
        return new Genotype(sample, Arrays.asList(a1, a2));
    }

    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError, double... pls) {
        return new Genotype(sample, Arrays.asList(a1, a2), log10pError, pls);
    }


    private Genotype makeG(String sample, Allele a1, Allele a2, double log10pError) {
        return new Genotype(sample, Arrays.asList(a1, a2), log10pError);
    }

    private VariantContext makeVC(String source, List<Allele> alleles) {
        return makeVC(source, alleles, null, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Genotype... g1) {
        return makeVC(source, alleles, Arrays.asList(g1));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, String filter) {
        return makeVC(source, alleles, filter.equals(".") ? null : new HashSet<String>(Arrays.asList(filter)));
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Set<String> filters) {
        return makeVC(source, alleles, null, filters);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes) {
        return makeVC(source, alleles, genotypes, null);
    }

    private VariantContext makeVC(String source, List<Allele> alleles, Collection<Genotype> genotypes, Set<String> filters) {
        int start = 10;
        int stop = start; // alleles.contains(ATC) ? start + 3 : start;
        return new VariantContext(source, "1", start, stop, alleles,
                genotypes == null ? null : VariantContext.genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes),
                1.0, filters, null, Cref.getBases()[0]);
    }

    // --------------------------------------------------------------------------------
    //
    // Test allele merging
    //
    // --------------------------------------------------------------------------------

    private class MergeAllelesTest extends TestDataProvider {
        List<List<Allele>> inputs;
        List<Allele> expected;

        private MergeAllelesTest(List<Allele>... arg) {
            super(MergeAllelesTest.class);
            LinkedList<List<Allele>> all = new LinkedList<List<Allele>>(Arrays.asList(arg));
            expected = all.pollLast();
            inputs = all;
        }

        public String toString() {
            return String.format("MergeAllelesTest input=%s expected=%s", inputs, expected);
        }
    }
    @DataProvider(name = "mergeAlleles")
    public Object[][] mergeAllelesData() {
        // first, do no harm
        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref));

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref),
                Arrays.asList(Aref));

        new MergeAllelesTest(Arrays.asList(Aref),
                Arrays.asList(Aref, T),
                Arrays.asList(Aref, T));

        new MergeAllelesTest(Arrays.asList(Aref, C),
                Arrays.asList(Aref, T),
                Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        new MergeAllelesTest(Arrays.asList(Aref, C, T),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, C, T), Arrays.asList(Aref, C, T));
        new MergeAllelesTest(Arrays.asList(Aref, T, C), Arrays.asList(Aref, T, C));

        new MergeAllelesTest(Arrays.asList(Aref, T, C),
                Arrays.asList(Aref, C),
                Arrays.asList(Aref, T, C)); // in order of appearence

        // The following is actually a pathological case - there's no way on a vcf to represent a null allele that's non-variant.
        // The code converts this (correctly) to a single-base non-variant vc with whatever base was there as a reference.
        new MergeAllelesTest(Arrays.asList(delRef),
                Arrays.asList(Cref));

        new MergeAllelesTest(Arrays.asList(delRef),
                Arrays.asList(delRef, ATC),
                Arrays.asList(delRef, ATC));

        new MergeAllelesTest(Arrays.asList(delRef),
                Arrays.asList(delRef, ATC, ATCATC),
                Arrays.asList(delRef, ATC, ATCATC));

        // alleles in the order we see them
        new MergeAllelesTest(Arrays.asList(delRef, ATCATC),
                Arrays.asList(delRef, ATC, ATCATC),
                Arrays.asList(delRef, ATCATC, ATC));

        // same
        new MergeAllelesTest(Arrays.asList(delRef, ATC),
                Arrays.asList(delRef, ATCATC),
                Arrays.asList(delRef, ATC, ATCATC));

        return MergeAllelesTest.getTests(MergeAllelesTest.class);
    }

    @Test(dataProvider = "mergeAlleles")
    public void testMergeAlleles(MergeAllelesTest cfg) {
        final List<VariantContext> inputs = new ArrayList<VariantContext>();

        int i = 0;
        for ( final List<Allele> alleles : cfg.inputs ) {
            final String name = "vcf" + ++i;
            inputs.add(makeVC(name, alleles));
        }

        final List<String> priority = vcs2priority(inputs);

        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                inputs, priority,
                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.PRIORITIZE, false, false, "set", false, false);

        Assert.assertEquals(merged.getAlleles(), cfg.expected);
    }

    // --------------------------------------------------------------------------------
    //
    // Test rsID merging
    //
    // --------------------------------------------------------------------------------

    private class SimpleMergeRSIDTest extends TestDataProvider {
        List<String> inputs;
        String expected;

        private SimpleMergeRSIDTest(String... arg) {
            super(SimpleMergeRSIDTest.class);
            LinkedList<String> allStrings = new LinkedList<String>(Arrays.asList(arg));
            expected = allStrings.pollLast();
            inputs = allStrings;
        }

        public String toString() {
            return String.format("SimpleMergeRSIDTest vc=%s expected=%s", inputs, expected);
        }
    }

    @DataProvider(name = "simplemergersiddata")
    public Object[][] createSimpleMergeRSIDData() {
        new SimpleMergeRSIDTest(".", ".");
        new SimpleMergeRSIDTest(".", ".", ".");
        new SimpleMergeRSIDTest("rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs1", "rs1");
        new SimpleMergeRSIDTest(".", "rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1,rs2");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1", "rs1,rs2"); // duplicates
        new SimpleMergeRSIDTest("rs2", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", "rs1", ".", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", ".", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs1", ".", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs3", "rs1,rs2,rs3");

        return SimpleMergeRSIDTest.getTests(SimpleMergeRSIDTest.class);
    }

    @Test(dataProvider = "simplemergersiddata")
    public void testRSIDMerge(SimpleMergeRSIDTest cfg) {
        final VariantContext snpVC1 = makeVC("snpvc1", Arrays.asList(Aref, T));
        final List<VariantContext> inputs = new ArrayList<VariantContext>();

        for ( final String id : cfg.inputs ) {
            MutableVariantContext vc = new MutableVariantContext(snpVC1);
            if ( ! id.equals(".") ) vc.setID(id);
            inputs.add(vc);
        }

        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                inputs, null,
                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.UNSORTED, false, false, "set", false, false);
        Assert.assertEquals(merged.getID(), cfg.expected.equals(".") ? null : cfg.expected);
    }

    // --------------------------------------------------------------------------------
    //
    // Test filtered merging
    //
    // --------------------------------------------------------------------------------

    private class MergeFilteredTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        String setExpected;
        VariantContextUtils.FilteredRecordMergeType type;


        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, String setExpected) {
            this(name, input1, input2, expected, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED, setExpected);
        }

        private MergeFilteredTest(String name, VariantContext input1, VariantContext input2, VariantContext expected, VariantContextUtils.FilteredRecordMergeType type, String setExpected) {
            super(MergeFilteredTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<VariantContext>(Arrays.asList(input1, input2));
            this.expected = expected;
            this.type = type;
            inputs = all;
            this.setExpected = setExpected;
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeFiltered")
    public Object[][] mergeFilteredData() {
        new MergeFilteredTest("AllPass",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                VariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("noFilters",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                VariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("oneFiltered",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", VariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("onePassOneFail",
                makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                String.format("1-%s2", VariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("AllFiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                VariantContextUtils.MERGE_FILTER_IN_ALL);

        // test ALL vs. ANY
        new MergeFilteredTest("FailOneUnfiltered",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "."),
                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                String.format("%s1-2", VariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("OneFailAllUnfilteredArg",
                makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                makeVC("2", Arrays.asList(Aref, T), "."),
                makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ALL_UNFILTERED,
                String.format("%s1-2", VariantContextUtils.MERGE_FILTER_PREFIX));

        // test excluding allele in filtered record
        new MergeFilteredTest("DontIncludeAlleleOfFilteredRecords",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                makeVC("3", Arrays.asList(Aref, T), "."),
                String.format("1-%s2", VariantContextUtils.MERGE_FILTER_PREFIX));

        // promotion of site from unfiltered to PASSES
        new MergeFilteredTest("UnfilteredPlusPassIsPass",
                makeVC("1", Arrays.asList(Aref, T), "."),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                VariantContextUtils.MERGE_INTERSECTION);

        new MergeFilteredTest("RefInAll",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                VariantContextUtils.MERGE_REF_IN_ALL);

        new MergeFilteredTest("RefInOne",
                makeVC("1", Arrays.asList(Aref), VariantContext.PASSES_FILTERS),
                makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                makeVC("3", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS),
                "2");

        return MergeFilteredTest.getTests(MergeFilteredTest.class);
    }

    @Test(dataProvider = "mergeFiltered")
    public void testMergeFiltered(MergeFilteredTest cfg) {
        final List<String> priority = vcs2priority(cfg.inputs);
        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                cfg.inputs, priority, cfg.type, VariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test set field
        Assert.assertEquals(merged.getAttribute("set"), cfg.setExpected);

        // test filter field
        Assert.assertEquals(merged.getFilters(), cfg.expected.getFilters());
    }

    // --------------------------------------------------------------------------------
    //
    // Test genotype merging
    //
    // --------------------------------------------------------------------------------

    private class MergeGenotypesTest extends TestDataProvider {
        List<VariantContext> inputs;
        VariantContext expected;
        List<String> priority;

        private MergeGenotypesTest(String name, String priority, VariantContext... arg) {
            super(MergeGenotypesTest.class, name);
            LinkedList<VariantContext> all = new LinkedList<VariantContext>(Arrays.asList(arg));
            this.expected = all.pollLast();
            inputs = all;
            this.priority = Arrays.asList(priority.split(","));
        }

        public String toString() {
            return String.format("%s input=%s expected=%s", super.toString(), inputs, expected);
        }
    }

    @DataProvider(name = "mergeGenotypes")
    public Object[][] mergeGenotypesData() {
        new MergeGenotypesTest("TakeGenotypeByPriority-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-1,2-nocall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, 1)));

        new MergeGenotypesTest("TakeGenotypeByPriority-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2)));

        new MergeGenotypesTest("NonOverlappingGenotypes", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, 2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1), makeG("s2", Aref, T, 2)));

        new MergeGenotypesTest("PreserveNoCall", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s2", Aref, T, 2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Allele.NO_CALL, Allele.NO_CALL, 1), makeG("s2", Aref, T, 2)));

        new MergeGenotypesTest("PerserveAlleles", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, C), makeG("s2", Aref, C, 2)),
                makeVC("3", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, 1), makeG("s2", Aref, C, 2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2), makeG("s3", Aref, T, 3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1), makeG("s3", Aref, T, 3)));

        new MergeGenotypesTest("TakeGenotypePartialOverlap-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2), makeG("s3", Aref, T, 3)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2), makeG("s3", Aref, T, 3)));

        //
        // merging genothpes with PLs
        //

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs", "1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1, 1, 2, 3)),
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1, 1, 2, 3)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles", "1",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6)));

        // first, do no harm
        new MergeGenotypesTest("OrderedPLs-3Alleles-2", "1",
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s2", Aref, C, 1, 1, 2, 3, 4, 5, 6)),
                makeVC("1", Arrays.asList(Aref, T, C), makeG("s1", Aref, T, 1, 1, 2, 3, 4, 5, 6), makeG("s2", Aref, C, 1, 1, 2, 3, 4, 5, 6)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-2,1", "2,1",
                makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2,4,0,2), makeG("s3", Aref, T, 3,3,0,2)),
                makeVC("3", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2,4,0,2), makeG("s3", Aref, T, 3,3,0,2)));

        new MergeGenotypesTest("TakeGenotypePartialOverlapWithPLs-1,2", "1,2",
                makeVC("1", Arrays.asList(Aref,ATC), makeG("s1", Aref, ATC, 1,5,0,3)),
                makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2,4,0,2), makeG("s3", Aref, T, 3,3,0,2)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, ATC, T), makeG("s1", Aref, ATC, 1), makeG("s3", Aref, T, 3)));

        new MergeGenotypesTest("MultipleSamplePLsDifferentOrder", "1,2",
                makeVC("1", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, 1, 1, 2, 3, 4, 5, 6)),
                makeVC("2", Arrays.asList(Aref, T, C), makeG("s2", Aref, T, 2, 6, 5, 4, 3, 2, 1)),
                // no likelihoods on result since type changes to mixed multiallelic
                makeVC("3", Arrays.asList(Aref, C, T), makeG("s1", Aref, C, 1), makeG("s2", Aref, T, 2)));

        return MergeGenotypesTest.getTests(MergeGenotypesTest.class);
    }

    @Test(dataProvider = "mergeGenotypes")
    public void testMergeGenotypes(MergeGenotypesTest cfg) {
        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                cfg.inputs, cfg.priority, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test genotypes
        assertGenotypesAreMostlyEqual(merged.getGenotypes(), cfg.expected.getGenotypes());
    }

    // necessary to not overload equals for genotypes
    private void assertGenotypesAreMostlyEqual(Map<String, Genotype> actual, Map<String, Genotype> expected) {
        if (actual == expected) {
            return;
        }

        if (actual == null || expected == null) {
            Assert.fail("Maps not equal: expected: " + expected + " and actual: " + actual);
        }

        if (actual.size() != expected.size()) {
            Assert.fail("Maps do not have the same size:" + actual.size() + " != " + expected.size());
        }

        for (Map.Entry<String, Genotype> entry : actual.entrySet()) {
            String key = entry.getKey();
            Genotype value = entry.getValue();
            Genotype expectedValue = expected.get(key);

            Assert.assertEquals(value.alleles, expectedValue.alleles, "Alleles in Genotype aren't equal");
            Assert.assertEquals(value.getNegLog10PError(), expectedValue.getNegLog10PError(), "GQ values aren't equal");
            Assert.assertEquals(value.hasLikelihoods(), expectedValue.hasLikelihoods(), "Either both have likelihoods or both not");
            if ( value.hasLikelihoods() )
                Assert.assertEquals(value.getLikelihoods().getAsVector(), expectedValue.getLikelihoods().getAsVector(), "Genotype likelihoods aren't equal");
        }
    }

    @Test
    public void testMergeGenotypesUniquify() {
        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1));
        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2));

        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                Arrays.asList(vc1, vc2), null, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.UNIQUIFY, false, false, "set", false, false);

        // test genotypes
        Assert.assertEquals(merged.getGenotypes().keySet(), new HashSet<String>(Arrays.asList("s1.1", "s1.2")));
    }

    @Test(expectedExceptions = UserException.class)
    public void testMergeGenotypesRequireUnique() {
        final VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), makeG("s1", Aref, T, 1));
        final VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), makeG("s1", Aref, T, 2));

        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                Arrays.asList(vc1, vc2), null, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE, false, false, "set", false, false);
    }

    // --------------------------------------------------------------------------------
    //
    // Misc. tests
    //
    // --------------------------------------------------------------------------------

    @Test
    public void testAnnotationSet() {
        for ( final boolean annotate : Arrays.asList(true, false)) {
            for ( final String set : Arrays.asList("set", "combine", "x")) {
                final List<String> priority = Arrays.asList("1", "2");
                VariantContext vc1 = makeVC("1", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS);
                VariantContext vc2 = makeVC("2", Arrays.asList(Aref, T), VariantContext.PASSES_FILTERS);

                final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                        Arrays.asList(vc1, vc2), priority, VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                        VariantContextUtils.GenotypeMergeType.PRIORITIZE, annotate, false, set, false, false);

                if ( annotate )
                    Assert.assertEquals(merged.getAttribute(set), VariantContextUtils.MERGE_INTERSECTION);
                else
                    Assert.assertFalse(merged.hasAttribute(set));
            }
        }
    }

    private static final List<String> vcs2priority(final Collection<VariantContext> vcs) {
        final List<String> priority = new ArrayList<String>();

        for ( final VariantContext vc : vcs ) {
            priority.add(vc.getSource());
        }

        return priority;
    }
}
