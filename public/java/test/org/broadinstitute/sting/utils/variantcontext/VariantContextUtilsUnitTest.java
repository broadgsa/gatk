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
    Allele Aref, T, delRef, ATC;
    Genotype ref1, snp1, snp2, indel1, indelref;
    private GenomeLocParser genomeLocParser;
    VariantContext refVC, snpVC1, snpVC2, snpVC3, snpVC4, indelVC1, indelVC2, indelVC3;

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
        ATC = Allele.create("ATC");

        ref1 = new Genotype("ref1", Arrays.asList(Aref, Aref), 5, new double[]{0, 5, 10});
        snp1 = new Genotype("snp1", Arrays.asList(Aref,T), 10, new double[]{10, 0, 20});
        snp2 = new Genotype("snp2", Arrays.asList(T,T), 15, new double[]{25, 15, 0});
        indelref = new Genotype("indelref", Arrays.asList(delRef,delRef), 25, new double[]{0, 25, 30});
        indel1 = new Genotype("indel1", Arrays.asList(delRef,ATC), 20, new double[]{20, 0, 30});

        refVC = makeVC("refvc", Arrays.asList(Aref), Arrays.asList(ref1));
        snpVC1 = makeVC("snpvc1", Arrays.asList(Aref, T), Arrays.asList(snp1));
        snpVC2 = makeVC("snpvc2", Arrays.asList(Aref, T), Arrays.asList(snp1, snp2));
        snpVC3 = makeVC("snpvc3", Arrays.asList(Aref, T), Arrays.asList(ref1, snp1));
        snpVC4 = makeVC("snpvc4", Arrays.asList(Aref, T), Arrays.asList(ref1, snp1, snp2));
        indelVC1 = makeVC("indelvc1", Arrays.asList(delRef), Arrays.asList(indelref));
        indelVC2 = makeVC("indelvc2", Arrays.asList(delRef, ATC), Arrays.asList(indel1));
        indelVC3 = makeVC("indelvc3", Arrays.asList(delRef, ATC), Arrays.asList(indelref, indel1));
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
                VariantContext.genotypeCollectionToMap(new TreeMap<String, Genotype>(), genotypes),
                1.0, filters, null, (byte)'C');
    }

    private class SimpleMergeTest extends TestDataProvider {
        List<VariantContext> inputVCs;
        VariantContext expectedVC;

        private SimpleMergeTest(VariantContext... vcsArg) {
            super(SimpleMergeTest.class);
            LinkedList<VariantContext> allVCs = new LinkedList<VariantContext>(Arrays.asList(vcsArg));
            expectedVC = allVCs.pollLast();
            inputVCs = allVCs;
        }

        public String toString() {
            return String.format("SimpleMergeTest vc=%s expected=%s", inputVCs, expectedVC);
        }
    }

    @DataProvider(name = "simplemergedata")
    public Object[][] createSimpleMergeData() {
        // first, do no harm
        new SimpleMergeTest(refVC,    refVC);
        new SimpleMergeTest(snpVC1,   snpVC1);
        new SimpleMergeTest(indelVC1, indelVC1);
        new SimpleMergeTest(indelVC3, indelVC3);

        new SimpleMergeTest(refVC,  snpVC1, snpVC3);
        new SimpleMergeTest(snpVC1, snpVC2, snpVC2);
        new SimpleMergeTest(refVC,  snpVC2, snpVC4);

        new SimpleMergeTest(indelVC1,  indelVC2, indelVC3);
        new SimpleMergeTest(indelVC1,  indelVC3, indelVC3);
        new SimpleMergeTest(indelVC2,  indelVC3, indelVC3);

        return SimpleMergeTest.getTests(SimpleMergeTest.class);
    }

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
        new SimpleMergeRSIDTest("rs1", "rs1");
        new SimpleMergeRSIDTest(".", "rs1", "rs1");
        new SimpleMergeRSIDTest("rs1", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs1,rs2");
        new SimpleMergeRSIDTest("rs2", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", "rs1", ".", "rs2,rs1");
        new SimpleMergeRSIDTest("rs2", ".", "rs1", "rs2,rs1");
        new SimpleMergeRSIDTest("rs1", ".", ".", "rs1");
        new SimpleMergeRSIDTest("rs1", "rs2", "rs3", "rs1,rs2,rs3");

        return SimpleMergeRSIDTest.getTests(SimpleMergeRSIDTest.class);
    }

    @Test(dataProvider = "simplemergersiddata")
    public void testRSIDMerge(SimpleMergeRSIDTest cfg) {
        List<VariantContext> inputs = new ArrayList<VariantContext>();
        for ( String id : cfg.inputs ) {
            MutableVariantContext vc = new MutableVariantContext(snpVC1);
            if ( ! id.equals(".") ) vc.setID(id);
            inputs.add(vc);

        }

        VariantContext merged = myMerge(inputs);
        Assert.assertEquals(merged.getID(), cfg.expected.equals(".") ? null : cfg.expected);
    }

    private VariantContext myMerge(List<VariantContext> inputs) {
        List<String> priority = new ArrayList<String>();
        for ( VariantContext vc : inputs ) priority.add(vc.getSource());

        return VariantContextUtils.simpleMerge(genomeLocParser,
                inputs, priority,
                VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED,
                VariantContextUtils.GenotypeMergeType.PRIORITIZE, true, false, "set", false, false);
    }

    // todo -- add tests for subset merging, especially with correct PLs
    // todo -- test priority list
    // todo -- test FilteredRecordMergeType
    // todo -- no annotate origin
    // todo -- test set key
    // todo -- test filtered are uncalled
}
