package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class VariantFiltrationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T VariantFiltration -o %s -NO_HEADER -R " + b36KGReference;
    }


    @Test
    public void testNoAction() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("b7b7c218e219cd923ce5b6eefc5b7171"));
        executeTest("test no action", spec);
    }

    @Test
    public void testClusteredSnps() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -window 10 --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6d45a19e4066e7de6ff6a61f43ffad2b"));
        executeTest("test clustered SNPs", spec);
    }

    @Test
    public void testMask1() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF3 " + validationDataLocation + "vcfexample2.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("65b5006bf3ee9d9d08a36d6b854773f2"));
        executeTest("test mask all", spec1);
    }

    @Test
    public void testMask2() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -maskName foo --mask:VCF " + validationDataLocation + "vcfMask.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("a275d36baca81a1ce03dbb528e95a069"));
        executeTest("test mask some", spec2);
    }

    @Test
    public void testMask3() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec3 = new WalkerTestSpec(
                baseTestString() + " -maskName foo -maskExtend 10 --mask:VCF " + validationDataLocation + "vcfMask.vcf --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("c9489e1c1342817c36ab4f0770609bdb"));
        executeTest("test mask extend", spec3);
    }

    @Test
    public void testFilter1() {
        WalkerTestSpec spec = new WalkerTestSpec(
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
                baseTestString() + " -filter 'DoC < 20 || FisherStrand > 20.0' -filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("327a611bf82c6c4ae77fbb6d06359f9d"));
        executeTest("test filter #1", spec);
    }

    @Test
    public void testFilter2() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -filter 'AlleleBalance < 70.0 && FisherStrand == 1.4' -filterName bar --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("7612b3460575402ad78fa4173178bdcc"));
        executeTest("test filter #2", spec);
    }

    @Test
    public void testFilterWithSeparateNames() {
        // note that this input if slightly malformed, but with the new properly
        // only when really needed genotype loading of VCF files we don't actually
        // fix the file in the output
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterName ABF -filter 'AlleleBalance < 0.7' --filterName FSF -filter 'FisherStrand == 1.4' --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("dce33441f58b284ac9ab94f8e64b84e3"));
        executeTest("test filter with separate names #2", spec);
    }

    @Test
    public void testGenotypeFilters1() {
        WalkerTestSpec spec1 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'GQ == 0.60' -G_filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("96b61e4543a73fe725e433f007260039"));
        executeTest("test genotype filter #1", spec1);
    }

    @Test
    public void testGenotypeFilters2() {
        WalkerTestSpec spec2 = new WalkerTestSpec(
                baseTestString() + " -G_filter 'AF == 0.04 && isHomVar == 1' -G_filterName foo --variant:VCF3 " + validationDataLocation + "vcfexample2.vcf -L 1:10,020,000-10,021,000", 1,
                Arrays.asList("6c8112ab17ce39c8022c891ae73bf38e"));
        executeTest("test genotype filter #2", spec2);
    }

    @Test
    public void testDeletions() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " --filterExpression 'QUAL < 100' --filterName foo --variant:VCF " + validationDataLocation + "twoDeletions.vcf", 1,
                Arrays.asList("569546fd798afa0e65c5b61b440d07ac"));
        executeTest("test deletions", spec);
    }
}
