package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    @Test
    public void testDiscordanceNoSampleSpecified() {
        String testFile = testDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant " + b37hapmapGenotypes + " -disc " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("133fd0ded0bb213097cbe68995afbb7e")
        );
        spec.disableShadowBCF();

        executeTest("testDiscordanceNoSampleSpecified--" + testFile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = testDir + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variant " + testfile),
                1,
                Arrays.asList("6c1a9e64a00a5b312531729bc73b5183")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = testDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant " + b37hapmapGenotypes + " -disc " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("f64c90c4cca470f1095d9fa2062eac3e")
        );
        spec.disableShadowBCF();

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile),
            1,
            Arrays.asList("eb1d0ff1db27413c14ea1af52b2f74c8")
        );
        spec.disableShadowBCF();
        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testSampleExclusion() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_sn A -xl_sf " + samplesFile + " --variant " + testfile,
            1,
            Arrays.asList("ed0f40334a82aa8e4698d5bfd8ed4d52")
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusion--" + testfile, spec);
    }


    @Test
    public void testConcordance() {
        String testFile = testDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc " + b37hapmapGenotypes + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("9da5dab3d344c1c0a5987b15e60fa082")
        );
        spec.disableShadowBCF();

        executeTest("testConcordance--" + testFile, spec);
    }

    @Test
    public void testVariantTypeSelection() {
        String testFile = testDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -restrictAllelesTo MULTIALLELIC -selectType MIXED --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("30b89b3a6706f7f46b23bfb3be69cc8e")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = testDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("4eeab0dd18712d0ca5ebe77c24ec989f")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRegenotype() {
        String testFile = testDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -regenotype -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("5bf9663274ceb552f5469f8c1dfc22ed")
        );

        executeTest("testRegenotype--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = testDir + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("cb9932f9a7aa2e53af605b30d88ad43f")
        );

        executeTest("testMultipleRecordsAtOnePosition--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = testDir + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("920605cc2182026e3f54c009f6a04141")
        );

        executeTest("testNoGTs--" + testFile, spec);
    }

    @Test(enabled = false)
    public void testParallelization2() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";
        WalkerTestSpec spec;

        spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile + " -nt 2"),
            1,
            Arrays.asList("eb1d0ff1db27413c14ea1af52b2f74c8")
        );
        spec.disableShadowBCF();
        executeTest("testParallelization (2 threads)--" + testfile, spec);
    }

    @Test(enabled = false)
    public void testParallelization4() {
            String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
            String samplesFile = validationDataLocation + "SelectVariants.samples.txt";
            WalkerTestSpec spec;
            spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile + " -nt 4"),
            1,
            Arrays.asList("eb1d0ff1db27413c14ea1af52b2f74c8")
        );
        spec.disableShadowBCF();

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }

    @Test
    public void testSelectFromMultiAllelic() {
        String testfile = testDir + "multi-allelic.bi-allelicInGIH.vcf";
        String samplesFile = testDir + "GIH.samples.list";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " -o %s --no_cmdline_in_header -sf " + samplesFile + " --excludeNonVariants --variant " + testfile,
                1,
                Arrays.asList("e3d2e00dc7bfff85b87f78b6162ad333")
        );
        executeTest("test select from multi allelic with excludeNonVariants --" + testfile, spec);
    }
}
