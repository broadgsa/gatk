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
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variants:VCF3 " + testfile + " -NO_HEADER"),
            1,
            Arrays.asList("d18516c1963802e92cb9e425c0b75fd6")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = validationDataLocation + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variants:VCF3 " + testfile + " -NO_HEADER"),
                1,
                Arrays.asList("b74038779fe6485dbb8734ae48178356")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variants:VCF " + b37hapmapGenotypes + " -disc:VCF " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("78e6842325f1f1bc9ab30d5e7737ee6e")
        );

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testConcordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc:VCF " + b37hapmapGenotypes + " --variants:VCF " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("d2ba3ea30a810f6f0fbfb1b643292b6a")
        );

        executeTest("testConcordance--" + testFile, spec);
    }

}
