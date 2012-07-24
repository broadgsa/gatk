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
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("d88bdae45ae0e74e8d8fd196627e612c")
        );
        spec.disableShadowBCF();

        executeTest("testDiscordanceNoSampleSpecified--" + testFile, spec);
    }

    @Test
    public void testRepeatedLineSelection() {
        String testfile = privateTestDir + "test.dup.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -sn A -sn B -sn C --variant " + testfile),
                1,
                Arrays.asList("3d98a024bf3aecbd282843e0af89d0e6")
        );

        executeTest("testRepeatedLineSelection--" + testfile, spec);
    }

    @Test
    public void testDiscordance() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 --variant "
                        + b37hapmapGenotypes + " -disc " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("54289033d35d32b8ebbb38c51fbb614c")
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
                Arrays.asList("4386fbb258dcef4437495a37f5a83c53")
        );
        spec.disableShadowBCF();
        executeTest("testComplexSelection--" + testfile, spec);
    }

    @Test
    public void testNonExistingFieldSelection() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(" -env -ef -select 'foo!=0||DP>0' --variant " + testfile),
                1,
                Arrays.asList("44e77cea624cfff2b8acc3a4b30485cb")    // should yield empty vcf because the foo!=0 will yield complete expression false
        );
        spec.disableShadowBCF();
        executeTest("testNonExistingSelection--" + testfile, spec);
    }

    @Test
    public void testSampleExclusion() {
        String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
        String samplesFile = validationDataLocation + "SelectVariants.samples.txt";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -L 1:1-1000000 -o %s --no_cmdline_in_header -xl_sn A -xl_sf " + samplesFile + " --variant " + testfile,
                1,
                Arrays.asList("1f5c72951a35667c4bdf1be153787e27")
        );
        spec.disableShadowBCF();

        executeTest("testSampleExclusion--" + testfile, spec);
    }


    @Test
    public void testConcordance() {
        String testFile = privateTestDir + "NA12878.hg19.example1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + hg19Reference + " -sn NA12878 -L 20:1012700-1020000 -conc "
                        + b37hapmapGenotypes + " --variant " + testFile
                        + " -o %s --no_cmdline_in_header -U LENIENT_VCF_PROCESSING",
                1,
                Arrays.asList("946e7f2e0ae08dc0e931c1634360fc46")
        );
        spec.disableShadowBCF();

        executeTest("testConcordance--" + testFile, spec);
    }

    @Test
    public void testVariantTypeSelection() {
        String testFile = privateTestDir + "complexExample1.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -restrictAllelesTo MULTIALLELIC -selectType MIXED --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ca2b70e3171420b08b0a2659bfe2a794")
        );

        executeTest("testVariantTypeSelection--" + testFile, spec);
    }

    @Test
    public void testUsingDbsnpName() {
        String testFile = privateTestDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("d12ae1617deb38f5ed712dc326935b9a")
        );

        executeTest("testUsingDbsnpName--" + testFile, spec);
    }

    @Test
    public void testRegenotype() {
        String testFile = privateTestDir + "combine.3.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -regenotype -sn NA12892 --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("c22ad8864d9951403672a24c20d6c3c2")
        );

        executeTest("testRegenotype--" + testFile, spec);
    }

    @Test
    public void testMultipleRecordsAtOnePosition() {
        String testFile = privateTestDir + "selectVariants.onePosition.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b36KGReference + " -select 'KG_FREQ < 0.5' --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("44f7c47395ca5b2afef5313f592c8cea")
        );

        executeTest("testMultipleRecordsAtOnePosition--" + testFile, spec);
    }

    @Test
    public void testNoGTs() {
        String testFile = privateTestDir + "vcf4.1.example.vcf";

        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " --variant " + testFile + " -o %s --no_cmdline_in_header",
                1,
                Arrays.asList("ef3c5f75074a5dd2b2cd2715856a2542")
        );

        executeTest("testNoGTs--" + testFile, spec);
    }

    @Test
    public void testSelectFromMultiAllelic() {
        String testfile = privateTestDir + "multi-allelic.bi-allelicInGIH.vcf";
        String samplesFile = privateTestDir + "GIH.samples.list";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SelectVariants -R " + b37KGReference + " -o %s --no_cmdline_in_header -sf " + samplesFile + " --excludeNonVariants --variant " + testfile,
                1,
                Arrays.asList("3ab35d5e81a29fb5db3e2add11c7e823")
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
        final String testFile = privateTestDir + "missingHeaderLine.vcf";
        final String cmd = "-T SelectVariants -R " + b36KGReference + " -sn NA12892 --variant:dbsnp "
                + testFile + " -o %s --no_cmdline_in_header"
                + (expectedException == null ? " -U LENIENT_VCF_PROCESSING" : "");
        WalkerTestSpec spec =
                expectedException != null
                        ? new WalkerTestSpec(cmd, 1, expectedException)
                        : new WalkerTestSpec(cmd, 1, Arrays.asList(""));
        spec.disableShadowBCF();

        executeTest(name, spec);
    }
}
