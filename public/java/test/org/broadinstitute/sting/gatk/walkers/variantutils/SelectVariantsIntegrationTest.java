package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectVariantsIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T SelectVariants -R " + b36KGReference + " -L 1 -o %s --no_cmdline_in_header" + args;
    }

    @Test
    public void testComplexSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile),
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
            "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_sn A -xl_sf " + samplesFile + " --variant " + testfile,
            1,
            Arrays.asList("730f021fd6ecf1d195dabbee2e233bfd")
        );

        executeTest("testSampleExclusion--" + testfile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = validationDataLocation + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variant " + testfile),
                1,
                Arrays.asList("5085a2f8cddfeae9f6274f905025184f")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant " + b37hapmapGenotypes + " -disc " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("929bbb96381541c162dc7e5462e26ea2")
        );

        executeTest("testDiscordance--" + testFile, spec);
    }

    @Test
    public void testDiscordanceNoSampleSpecified() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant " + b37hapmapGenotypes + " -disc " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("5d7d899c0c4954ec59104aebfe4addd5")
        );

        executeTest("testDiscordanceNoSampleSpecified--" + testFile, spec);
    }

    @Test
    public void testConcordance() {
        String testFile = validationDataLocation + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc " + b37hapmapGenotypes + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("d2ba3ea30a810f6f0fbfb1b643292b6a")
        );

        executeTest("testConcordance--" + testFile, spec);
    }

    @Test
    public void testVariantTypeSelection() {
        String testFile = validationDataLocation + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -restrictAllelesTo MULTIALLELIC -selectType MIXED --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("e0b12c0b47a8a7a988b3587b47bfa8cf")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = validationDataLocation + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("167a1265df820978a74c267df44d5c43")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRegenotype() {
        String testFile = validationDataLocation + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -regenotype -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("0fd8e52bdcd1f4b921d8fb5c689f196a")
        );

        executeTest("testRegenotype--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = validationDataLocation + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ffa2524380d84a870d2e4a33d9f3d31a")
        );

        executeTest("testMultipleRecordsAtOnePositionFirstIsFiltered--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = validationDataLocation + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("f17885e5cbd5387edb99112047ea43c1")
        );

        executeTest("testMultipleRecordsAtOnePositionFirstIsFiltered--" + testFile, spec);
    }

    @Test
    public void testParallelization() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";
        WalkerTestSpec spec;

        spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile + " -nt 2"),
            1,
            Arrays.asList("d18516c1963802e92cb9e425c0b75fd6")
        );
        executeTest("testParallelization (2 threads)--" + testfile, spec);

        spec = new WalkerTestSpec(
            baseTestString(" -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile + " -nt 4"),
            1,
            Arrays.asList("d18516c1963802e92cb9e425c0b75fd6")
        );

        executeTest("testParallelization (4 threads)--" + testfile, spec);
    }

    @Test
    public void testSelectFromMultiAllelic() {
        String testfile = validationDataLocation + "multi-allelic.bi-allelicInGIH.vcf";
        String samplesFile = validationDataLocation + "GIH.samples.list";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " -o %s --no_cmdline_in_header -sf " + samplesFile + " --excludeNonVariants --variant " + testfile,
                1,
                Arrays.asList("3fb50cc1c955491048108956d7087c35")
        );
        executeTest("test select from multi allelic with excludeNonVariants --" + testfile, spec);
    }
}
