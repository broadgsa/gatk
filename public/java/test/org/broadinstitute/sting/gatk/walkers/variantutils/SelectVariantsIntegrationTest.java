package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
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
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header --allowMissingVCFHeaders --allowMissingVCFHeaders",
                1,
                Arrays.asList("18ca4d3e69874c5e55377ca2ad7c6014")
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
                Arrays.asList("ddb082f5f15944b7545b45c20440a474")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = testDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header --allowMissingVCFHeaders",
                1,
                Arrays.asList("2e0eaed94e2ab14349b967dc89a09823")
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
                Arrays.asList("88c142737080847e44dc8e8c982cfc74")
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
                Arrays.asList("338bfc635667793b89c349446982a36d")
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusion--" + testfile, spec);
    }


    @Test
    public void testConcordance() {
        String testFile = testDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc "
                        + b37hapmapGenotypes + " --variant " + testFile
                        + " -o %s --no_cmdline_in_header --allowMissingVCFHeaders",
                1,
                Arrays.asList("7d2c7deeacf307d140a85a34629c6539")
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
                Arrays.asList("f0509f4c2c6c9eff3bec8171f00762b3")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = testDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("9162a67ccb4201c0542f30d14967f2d5")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRegenotype() {
        String testFile = testDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -regenotype -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("be38bdc7bd88f5d09cf1a9d55cfecb0b")
        );

        executeTest("testRegenotype--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = testDir + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("7dabe6910a40a9db5006c6af66b4617c")
        );

        executeTest("testMultipleRecordsAtOnePosition--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = testDir + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("034965918283d5597bb443e65cd20cf8")
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
                Arrays.asList("88c142737080847e44dc8e8c982cfc74")
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
                Arrays.asList("88c142737080847e44dc8e8c982cfc74")
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
                Arrays.asList("3cbabcd256070ea113ddc04fef3b6900")
        );
        executeTest("test select from multi allelic with excludeNonVariants --" + testfile, spec);
    }

    @Test()
    public void testFileWithoutInfoLineInHeader() {
        testFileWithoutInfoLineInHeader("testFileWithoutInfoLineInHeader", UserException.class);
    }

    @Test()
    public void testFileWithoutInfoLineInHeaderWithOverride() {
        testFileWithoutInfoLineInHeader("testFileWithoutInfoLineInHeaderWithOverride", null);
    }

    private void testFileWithoutInfoLineInHeader(final String name, final Class expectedException) {
        final String testFile = testDir + "missingHeaderLine.vcf";
        final String cmd = "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp "
                + testFile + " -o %s --no_cmdline_in_header"
                + (expectedException == null ? " -allowMissingVCFHeaders" : "");
        WalkerTestSpec spec =
                expectedException != null
                        ? new WalkerTestSpec(cmd, 1, expectedException)
                        : new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();

        executeTest(name, spec);
    }
}
