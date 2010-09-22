package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

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
            baseTestString(" -sn A -sn '[CDH]' -sn " + samplesFile + " -env -ef -select 'AF < 0.2' -B:variant,VCF " + testfile + " -NO_HEADER"),
            1,
            Arrays.asList("3a15628b5980031c629c0c33e7e60b40")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }
}
