package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VCFIntegrationTest extends WalkerTest {

    @Test(enabled = true)
    public void testReadingAndWritingWitHNoChanges() {

        String md5ofInputVCF = "d991abe6c6a7a778a60a667717903be0";
        String testVCF = privateTestDir + "vcf4.1.example.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T VariantAnnotator --variant " + testVCF + " -L " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList(md5ofInputVCF));
        List<File> result = executeTest("Test Variant Annotator with no changes", spec1).getFirst();

        String test2 = baseCommand + "-T VariantsToVCF --variant " + result.get(0).getAbsolutePath();
        WalkerTestSpec spec2 = new WalkerTestSpec(test2, 1, Arrays.asList(md5ofInputVCF));
        executeTest("Test Variants To VCF from new output", spec2);
    }

    @Test(enabled = true)
    // See https://getsatisfaction.com/gsa/topics/support_vcf_4_1_structural_variation_breakend_alleles?utm_content=topic_link&utm_medium=email&utm_source=new_topic
    public void testReadingAndWritingBreakpointAlleles() {
        String testVCF = privateTestDir + "breakpoint-example.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("13329ba7360a8beb3afc02569e5a20c4"));
        executeTest("Test reading and writing breakpoint VCF", spec1);
    }

    @Test
    public void testReadingAndWritingSamtools() {
        String testVCF = privateTestDir + "samtools.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("0f82ac11852e7f958c1a0ce52398c2ae"));
        executeTest("Test reading and writing samtools vcf", spec1);
    }

    @Test
    public void testWritingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.vcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("9773d6a121cfcb18d090965bc520f120"));
        executeTest("Test writing samtools WEx BCF example", spec1);
    }

    @Test(enabled = false)
    public void testReadingSamtoolsWExBCFExample() {
        String testVCF = privateTestDir + "ex2.bcf";
        String baseCommand = "-R " + b36KGReference + " --no_cmdline_in_header -o %s ";
        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("63a2e0484ae37b0680514f53e0bf0c94"));
        executeTest("Test reading samtools WEx BCF example", spec1);
    }

    //
    //
    // Tests to ensure that -U LENIENT_VCF_PROCESS and header repairs are working
    //
    //

    @Test
    public void testFailingOnVCFWithoutHeaders() {
        runVCFWithoutHeaders("", "", UserException.class, false);
    }

    @Test
    public void testPassingOnVCFWithoutHeadersWithLenientProcessing() {
        runVCFWithoutHeaders("-U LENIENT_VCF_PROCESSING", "6de8cb7457154dd355aa55befb943f88", null, true);
    }

    @Test
    public void testPassingOnVCFWithoutHeadersRepairingHeaders() {
        runVCFWithoutHeaders("-repairVCFHeader " + privateTestDir + "vcfexample2.justHeader.vcf", "ff61e9cad6653c7f93d82d391f7ecdcb", null, false);
    }

    private void runVCFWithoutHeaders(final String moreArgs, final String expectedMD5, final Class expectedException, final boolean disableBCF) {
        final String testVCF = privateTestDir + "vcfexample2.noHeader.vcf";
        final String baseCommand = "-R " + b37KGReference
                + " --no_cmdline_in_header -o %s "
                + "-T VariantsToVCF -V " + testVCF + " " + moreArgs;
        WalkerTestSpec spec1 = expectedException != null
                ? new WalkerTestSpec(baseCommand, 1, expectedException)
                : new WalkerTestSpec(baseCommand, 1, Arrays.asList(expectedMD5));
        if ( disableBCF )
            spec1.disableShadowBCF();
        executeTest("Test reading VCF without header lines with additional args " + moreArgs, spec1);
    }
}
