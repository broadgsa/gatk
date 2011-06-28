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
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' -B:variant,VCF3 " + testfile + " -NO_HEADER"),
            1,
            Arrays.asList("1b9d551298dc048c7d36b60440ff4d50")
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

    @Test
    public void testDiscordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -disc myvar -L 20:1012700-1020000 -B:variant,VCF " + b37hapmapGenotypes + " -B:myvar,VCF " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("97621ae8f29955eedfc4e0be3515fcb9")
        );

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testConcordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -conc hapmap -L 20:1012700-1020000 -B:hapmap,VCF " + b37hapmapGenotypes + " -B:variant,VCF " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("a0ae016fdffcbe7bfb99fd3dbc311407")
        );

        executeTest("testConcordance--" + testFile, spec);
    }

}
