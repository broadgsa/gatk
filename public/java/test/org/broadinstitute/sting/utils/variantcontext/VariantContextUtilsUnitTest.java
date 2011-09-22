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

//    VariantContext refVC, snpVC1, snpVC2, snpVC3, snpVC4, indelVC1, indelVC2, indelVC3;
//    ref1 = new Genotype("ref1", Arrays.asList(Aref, Aref), 5, new double[]{0, 5, 10});
//    snp1 = new Genotype("snp1", Arrays.asList(Aref,T), 10, new double[]{10, 0, 20});
//    snp2 = new Genotype("snp2", Arrays.asList(T,T), 15, new double[]{25, 15, 0});
//    indelref = new Genotype("indelref", Arrays.asList(delRef,delRef), 25, new double[]{0, 25, 30});
//    indel1 = new Genotype("indel1", Arrays.asList(delRef,ATC), 20, new double[]{20, 0, 30});
//
//    refVC = makeVC("refvc", Arrays.asList(Aref), Arrays.asList(ref1));
//    snpVC1 = makeVC("snpvc1", Arrays.asList(Aref, T), Arrays.asList(snp1));
//    snpVC2 = makeVC("snpvc2", Arrays.asList(Aref, T), Arrays.asList(snp1, snp2));
//    snpVC3 = makeVC("snpvc3", Arrays.asList(Aref, T), Arrays.asList(ref1, snp1));
//    snpVC4 = makeVC("snpvc4", Arrays.asList(Aref, T), Arrays.asList(ref1, snp1, snp2));
//    indelVC1 = makeVC("indelvc1", Arrays.asList(delRef), Arrays.asList(indelref));
//    indelVC2 = makeVC("indelvc2", Arrays.asList(delRef, ATC), Arrays.asList(indel1));
//    indelVC3 = makeVC("indelvc3", Arrays.asList(delRef, ATC), Arrays.asList(indelref, indel1));

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

    // todo -- add tests for subset merging, especially with correct PLs
    // todo -- test priority list: intersection, filtered in all, reference in all, X-filteredInX, X
    // todo -- test FilteredRecordMergeType
    // todo -- no annotate origin
    // todo -- test set key
}
