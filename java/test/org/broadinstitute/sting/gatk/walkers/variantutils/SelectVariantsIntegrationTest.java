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
            baseTestString(" -sn A -sn '[CDH]' -sn " + samplesFile + " -env -ef -select 'DP < 250' -B:variant,VCF3 " + testfile + " -NO_HEADER"),
            1,
            Arrays.asList("1df84e2b755cce19f1876710ec38dd2c")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = validationDataLocation + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C -B:variant,VCF3 " + testfile + " -NO_HEADER"),
                1,
                Arrays.asList("5ba7536a0819421b330350a160e4261a")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }
}
