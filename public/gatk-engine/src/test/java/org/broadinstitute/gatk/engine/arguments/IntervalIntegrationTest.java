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

package org.broadinstitute.gatk.engine.arguments;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

/**
 * Test the GATK core interval parsing mechanism.
 */
public class IntervalIntegrationTest extends WalkerTest {
    @Test(enabled = true)
    public void testAllImplicitIntervalParsing() {
        String md5 = "7821db9e14d4f8e07029ff1959cd5a99";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testAllIntervalsImplicit",spec);
    }

    @Test
    public void testUnmappedReadInclusion() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestPrintReadsWalker" +
                        " -I " + validationDataLocation + "MV1994.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -L unmapped" +
                        " -U",
                        0, // two output files
                        Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadInclusion",".bam");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("c66bb2c3c5382e2acff09b2b359562bb",createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("fadcdf88597b9609c5f2a17f4c6eb455", createTempFileFromBase(baseOutputFile.getAbsolutePath().substring(0,baseOutputFile.getAbsolutePath().indexOf(".bam"))+".bai"));

        executeTest("testUnmappedReadInclusion",spec);
    }

    @Test
    public void testMultipleIntervalInclusionOnCRAM() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestPrintReadsWalker" +
                        " -I " + validationDataLocation + "MV1994.cram" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -L Escherichia_coli_K12:11000" +
                        " -L Escherichia_coli_K12:12000" +
                        " -L Escherichia_coli_K12:13000" +
                        " -L Escherichia_coli_K12:14000" +
                        " -L Escherichia_coli_K12:15000" +
                        " -L Escherichia_coli_K12:16000" +
                        " -L Escherichia_coli_K12:17000" +
                        " -L unmapped",
                0, // two output files
                Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadInclusion", ".cram");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("0f11cc035455cd68fb388e33aaf5feff", createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("ebbe6e311b6bb240554ec96ed9809216", createTempFileFromBase(baseOutputFile.getAbsolutePath() + ".bai"));

        executeTest("testUnmappedReadInclusionCRAM", spec);
    }

    @Test
    public void testMixedMappedAndUnmapped() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestPrintReadsWalker" +
                        " -I " + validationDataLocation + "MV1994.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -L Escherichia_coli_K12:4630000-4639675" +
                        " -L unmapped" +
                        " -U",
                        0, // two output files
                        Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadInclusion",".bam");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("c64cff3ed376bc8f2977078dbdac4518",createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("fa90ff91ac0cc689c71a3460a3530b8b", createTempFileFromBase(baseOutputFile.getAbsolutePath().substring(0,baseOutputFile.getAbsolutePath().indexOf(".bam"))+".bai"));

        executeTest("testUnmappedReadInclusion",spec);
    }


    @Test(enabled = false)
    public void testUnmappedReadExclusion() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestPrintReadsWalker" +
                        " -I " + validationDataLocation + "MV1994.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -XL unmapped" +
                        " -U",
                        0, // two output files
                        Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadExclusion",".bam");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("80887ba488e53dabd9596ff93070ae75",createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("b341d808ecc33217f37c0c0cde2a3e2f", createTempFileFromBase(baseOutputFile.getAbsolutePath().substring(0,baseOutputFile.getAbsolutePath().indexOf(".bam"))+".bai"));

        executeTest("testUnmappedReadExclusion",spec);
    }

    @Test(enabled = true)
    public void testIntervalParsingFromFile() {
        String md5 = "48a24b70a0b376535542b996af517398";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testIntervalParsingFromFile", spec);
    }

    @Test(enabled = true)
    public void testIntervalMergingFromFiles() {
        String md5 = "9ae0ea9e3c9c6e1b9b6252c8395efdc1";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf" +
                        " -L " + validationDataLocation + "intervalTest.2.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testIntervalMergingFromFiles", spec);
    }

    @Test(enabled = true)
    public void testIntervalExclusionsFromFiles() {
        String md5 = "26ab0db90d72e28ad0ba1e22ee510510";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf" +
                        " -XL " + validationDataLocation + "intervalTest.2.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testIntervalExclusionsFromFiles", spec);
    }

    @Test(enabled = true)
    public void testMixedIntervalMerging() {
        String md5 = "7c5aba41f53293b712fd86d08ed5b36e";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf" +
                        " -L chr1:1677524-1677528",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testMixedIntervalMerging", spec);
    }

    @Test(enabled = true)
    public void testBed() {
        String md5 = "cf4278314ef8e4b996e1b798d8eb92cf";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.bed",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testBed", spec);
    }

    @Test(enabled = true)
    public void testComplexVCF() {
        String md5 = "166d77ac1b46a1ec38aa35ab7e628ab5";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.3.vcf",
                1, // just one output file
                Collections.singletonList(md5));
        executeTest("testComplexVCF", spec);
    }

    @Test(enabled = true)
    public void testComplexVCFWithPadding() {
        String md5 = "649ee93d50739c656e94ec88a32c7ffe";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " --interval_padding 2" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.3.vcf",
                1, // just one output file
                Collections.singletonList(md5));
        executeTest("testComplexVCFWithPadding", spec);
    }

    @Test(enabled = true)
    public void testMergingWithComplexVCF() {
        String md5 = "6d7fce9fee471194aa8b5b6e47267f03";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf" +
                        " -XL " + validationDataLocation + "intervalTest.3.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testMergingWithComplexVCF", spec);
    }

    @Test(enabled = true)
    public void testEmptyVCF() {
        String md5 = "897316929176464ebc9ad085f31e7284";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.empty.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testEmptyVCFWarning", spec);
    }

    @Test(enabled = true)
    public void testIncludeExcludeIsTheSame() {
        String md5 = "897316929176464ebc9ad085f31e7284";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s" +
                        " -L " + validationDataLocation + "intervalTest.1.vcf" +
                        " -XL " + validationDataLocation + "intervalTest.1.vcf",
                        1, // just one output file
                        Collections.singletonList(md5));
        executeTest("testIncludeExcludeIsTheSame", spec);
    }

    @Test(enabled = true)
    public void testSymbolicAlleles() {
        String md5 = "52745056d2fd5904857bbd4984c08098";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T TestCountLociWalker" +
                        " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam" +
                        " -R " + b36KGReference +
                        " -o %s" +
                        " -L " + privateTestDir + "symbolic_alleles_1.vcf",
                1, // just one output file
                Collections.singletonList(md5));
        executeTest("testSymbolicAlleles", spec);
    }

    @Test
    public void testIntersectionOfLexicographicallySortedIntervals() {
        final String md5 = "18be9375e5a753f766616a51eb6131f0";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                " -T TestCountLociWalker" +
                " -I " + privateTestDir + "NA12878.4.snippet.bam" +
                " -R " + b37KGReference +
                " -L " + privateTestDir + "lexicographicallySortedIntervals.bed" +
                " -L 4" +
                " -isr INTERSECTION" +
                " -o %s",
                1, // just one output file
                Collections.singletonList(md5));
        executeTest("testIntersectionOfLexicographicallySortedIntervals", spec);
    }
}
