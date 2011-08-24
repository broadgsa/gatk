package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s -NO_HEADER" + args;
    }

    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant:VCF3 " + testfile),
            1,
            Arrays.asList("d18516c1963802e92cb9e425c0b75fd6")
        );

        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testSampleExclusion() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s -NO_HEADER -xl_sn A -xl_sf " + samplesFile + " --variant:VCF3 " + testfile,
            1,
            Arrays.asList("730f021fd6ecf1d195dabbee2e233bfd")
        );

        executeTest("testSampleExclusion--" + testfile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = validationDataLocation + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variant:VCF3 " + testfile),
                1,
                Arrays.asList("b74038779fe6485dbb8734ae48178356")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant:VCF " + b37hapmapGenotypes + " -disc:VCF " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("78e6842325f1f1bc9ab30d5e7737ee6e")
        );

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testConcordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc:VCF " + b37hapmapGenotypes + " --variant " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("d2ba3ea30a810f6f0fbfb1b643292b6a")
        );

        executeTest("testConcordance--" + testFile, spec);
    }

    @Test
    public void testVariantTypeSelection() {
        String testFile = validationDataLocation + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -restrictAllelesTo MULTIALLELIC -selectType MIXED --variant " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("e0b12c0b47a8a7a988b3587b47bfa8cf")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    @Test(enabled=false)
    public void testRemovePLs() {
        String testFile = validationDataLocation + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant " + testFile + " -o %s -NO_HEADER",
                1,
                Arrays.asList("")
        );

        executeTest("testWithPLs--" + testFile, spec);
    }
}
