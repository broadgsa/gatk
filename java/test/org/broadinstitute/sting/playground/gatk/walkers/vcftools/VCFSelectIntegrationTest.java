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
                Arrays.asList("ff15451c6138c75ece99c3bac91a9d4f"));
        executeTest("testVCFSelect1", spec);
    }

    @Test
    public void testVCFSelect2() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -L 1:10001290-10048590 ", 1,
                Arrays.asList("226ce18354d56f69d9506e7ae70e4eb9"));
        executeTest("testVCFSelect2", spec);
    }

    @Test
    public void testVCFSelectOr() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6' -match 'AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("62f6341ac676a919497784f792d3e22f"));
        executeTest("testVCFSelectOr", spec);
    }

    @Test
    public void testVCFSelectAnd() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B variant,VCF," + validationDataLocation + "vcfexample3.vcf -match 'HomopolymerRun == 6 && AF == 0.50' -L 1:10001290-10048590 ", 1,
                Arrays.asList("914845500749fbc1863d8226f31c96b3"));
        executeTest("testVCFSelectAnd", spec);
    }
}