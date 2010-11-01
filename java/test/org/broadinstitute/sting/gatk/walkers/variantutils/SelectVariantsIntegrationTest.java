package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s" + args;
    }

    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString(" -sn A -sn '[CDH]' -sn " + samplesFile + " -env -ef -select 'DP < 250' -B:variant,VCF " + testfile + " -NO_HEADER"),
            1,
            Arrays.asList("e22db73fd99fd4f28d472d407e40ead5")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = validationDataLocation + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C -B:variant,VCF " + testfile + " -NO_HEADER"),
                1,
                Arrays.asList("fae9822d5f7ad6c76b411e8ca0886409")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }
}
