package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VCFSelectIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VCFSelect -o %s -R " + oneKGLocation + "reference/human_b36_both.fasta";
    }


    @Test
    public void testVCFSelect1() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("b49ba344471444077bc6fe3c17e7bc3f"));
        executeTest("testVCFSelect1", spec);
    }

    @Test
    public void testVCFSelect2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -L 1:10001290-10048590 ", 1,
                Arrays.asList("517b4ae7058c3125ad6846c33a1a2e57"));
        executeTest("testVCFSelect2", spec);
    }

    @Test
    public void testVCFSelectOr() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("d77d8f938a61abd60fc813ff1a06bb0c"));
        executeTest("testVCFSelectOr", spec);
    }

    @Test
    public void testVCFSelectAnd() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6 && AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("ef05fc766482ffade95f1bbdb777770d"));
        executeTest("testVCFSelectAnd", spec);
    }
}