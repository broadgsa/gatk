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

// our package
package org.broadinstitute.sting.utils.variantcontext;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Priority;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.tracks.FeatureManager;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.yaml.snakeyaml.Yaml;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;


public class VariantContextUtilsUnitTest extends BaseTest {
    Allele Aref, T, C, delRef, ATC, ATCATC;
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
        delRef = Allele.create("-", true);
        T = Allele.create("T");
        C = Allele.create("C");
        ATC = Allele.create("ATC");
        ATCATC = Allele.create("ATCATC");
    }

    private VariantContext makeVC(String source, List<Allele> alleles) {
        return makeVC(source, alleles, null, null);
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
                1.0, filters, null, (byte)'C');
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
                             Arrays.asList(Aref, C, T)); // sorted by allele

        new MergeAllelesTest(Arrays.asList(Aref, C, T),
                             Arrays.asList(Aref, C),
                             Arrays.asList(Aref, C, T));

        new MergeAllelesTest(Arrays.asList(Aref, T, C),
                             Arrays.asList(Aref, C),
                             Arrays.asList(Aref, C, T)); // sorted by allele

        new MergeAllelesTest(Arrays.asList(delRef),
                             Arrays.asList(delRef)); // todo -- FIXME me GdA

        new MergeAllelesTest(Arrays.asList(delRef),
                             Arrays.asList(delRef, ATC),
                             Arrays.asList(delRef, ATC));

        new MergeAllelesTest(Arrays.asList(delRef),
                             Arrays.asList(delRef, ATC, ATCATC),
                             Arrays.asList(delRef, ATC, ATCATC));

        new MergeAllelesTest(Arrays.asList(delRef, ATCATC),
                             Arrays.asList(delRef, ATC, ATCATC),
                             Arrays.asList(delRef, ATC, ATCATC));

        new MergeAllelesTest(Arrays.asList(delRef, ATC),
                             Arrays.asList(delRef, ATCATC),
                             Arrays.asList(delRef, ATC, ATCATC));

        return MergeAllelesTest.getTests(MergeAllelesTest.class);
    }

    @Test(dataProvider = "mergeAlleles")
    public void testMergeAlleles(MergeAllelesTest cfg) {
        final List<String> priority = new ArrayList<String>();
        final List<VariantContext> inputs = new ArrayList<VariantContext>();

        int i = 0;
        for ( final List<Allele> alleles : cfg.inputs ) {
            final String name = "vcf" + ++i;
            priority.add(name);
            inputs.add(makeVC(name, alleles));
        }

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

        new MergeFilteredTest("FailOneUnfiltered",
                              makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                              makeVC("2", Arrays.asList(Aref, T), "."),
                              makeVC("3", Arrays.asList(Aref, T), "."),
                              String.format("%s1-2", VariantContextUtils.MERGE_FILTER_PREFIX));

        new MergeFilteredTest("AllFiltered",
                              makeVC("1", Arrays.asList(Aref, T), "FAIL"),
                              makeVC("2", Arrays.asList(Aref, T), "FAIL"),
                              makeVC("3", Arrays.asList(Aref, T), "FAIL"),
                              VariantContextUtils.MERGE_FILTER_IN_ALL);

        // test ALL vs. ANY
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

        return MergeFilteredTest.getTests(MergeFilteredTest.class);
    }

    @Test(dataProvider = "mergeFiltered")
    public void testMergeFiltered(MergeFilteredTest cfg) {
        final List<String> priority = new ArrayList<String>();

        for ( final VariantContext vc : cfg.inputs ) {
            priority.add(vc.getSource());
        }

        final VariantContext merged = VariantContextUtils.simpleMerge(genomeLocParser,
                cfg.inputs, priority, cfg.type, VariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);

        // test alleles are equal
        Assert.assertEquals(merged.getAlleles(), cfg.expected.getAlleles());

        // test set field
        Assert.assertEquals(merged.getAttribute("set"), cfg.setExpected);

        // test filter field
        Assert.assertEquals(merged.getFilters(), cfg.expected.getFilters());
    }


    // todo -- add tests for subset merging, especially with correct PLs
    // todo -- test priority list: intersection, filtered in all, reference in all, X-filteredInX, X
    // todo -- test FilteredRecordMergeType
    // todo -- no annotate origin
    // todo -- test set key
}
