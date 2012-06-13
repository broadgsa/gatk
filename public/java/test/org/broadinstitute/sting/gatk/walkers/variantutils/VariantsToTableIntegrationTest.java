/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.util.*;

public class VariantsToTableIntegrationTest extends WalkerTest {
    private String variantsToTableCmd(String moreArgs) {
        return "-R " + hg18Reference +
                " --variant:vcf " + testDir + "soap_gatk_annotated.vcf" +
                " -T VariantsToTable" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F TRANSITION -F DP -F SB -F set -F RankSumP -F refseq.functionalClass*" +
                " -L chr1 -o %s" + moreArgs;
    }

    private String variantsToTableMultiAllelicCmd(String moreArgs) {
        return "-R " + b37KGReference +
                " --variant " + testDir + "multiallelic.vcf" +
                " -T VariantsToTable" +
                " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                " -o %s" + moreArgs;
    }

    @Test(enabled = true)
    public void testComplexVariantsToTable() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableCmd(" -AMD"),
                Arrays.asList("e8f771995127b727fb433da91dd4ee98"));
        executeTest("testComplexVariantsToTable", spec);
    }

    @Test(enabled = true)
    public void testComplexVariantsToTableFail() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableCmd(""), 1, UserException.class);
        executeTest("testComplexVariantsToTable-FAIL", spec);
    }

    @Test(enabled = true)
    public void testMultiAllelicOneRecord() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableMultiAllelicCmd(""),
                Arrays.asList("13dd36c08be6c800f23988e6000d963e"));
        executeTest("testMultiAllelicOneRecord", spec);
    }

    @Test(enabled = true)
    public void testMultiAllelicSplitRecords() {
        WalkerTestSpec spec = new WalkerTestSpec(variantsToTableMultiAllelicCmd(" -SMA"),
                Arrays.asList("17a0fc80409d2fc00ad2bbb94b3a346b"));
        executeTest("testMultiAllelicSplitRecords", spec);
    }

    @Test(enabled = true)
    public void testGenotypeFields() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                " --variant " + testDir + "vcfexample2.vcf" +
                " -T VariantsToTable" +
                " -GF RD" +
                " -o %s",
                1,
                Arrays.asList("d43562e9b94f0e8e337d38a6829671ee"));
        executeTest("testGenotypeFields", spec);
    }

    @Test(enabled = true)
    public void testGenotypeFieldsWithInline() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant " + testDir + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -GF RD -GF GT -GF GQ" +
                        " -o %s",
                1,
                Arrays.asList("29744059742ae71fd6aabd29e5c391fb"));
        executeTest("testGenotypeFieldsWithInline", spec);
    }

    @Test(enabled = true)
    public void testMoltenOutput() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant " + testDir + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER" +
                        " --moltenize" +
                        " -o %s",
                1,
                Arrays.asList("30047a5e78a7f523bd2872ac8baccc0e"));
        executeTest("testMoltenOutput", spec);
    }

    @Test(enabled = true)
    public void testMoltenOutputWithGenotypeFields() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant " + testDir + "vcfexample2.vcf" +
                        " -T VariantsToTable" +
                        " -GF RD" +
                        " --moltenize" +
                        " -o %s",
                1,
                Arrays.asList("1d97fe63c249a995df4ce666382872d8"));
        executeTest("testMoltenOutputWithGenotypeFields", spec);
    }

    @Test(enabled = true)
    public void testMoltenOutputWithMultipleAlleles() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b37KGReference +
                        " --variant " + testDir + "multiallelic.vcf" +
                        " -T VariantsToTable" +
                        " -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F MULTI-ALLELIC -F AC -F AF" +
                        " --moltenize -SMA" +
                        " -o %s",
                1,
                Arrays.asList("c131e2c3cfb673c456cb160bda476101"));
        executeTest("testMoltenOutputWithMultipleAlleles", spec);
    }
}
