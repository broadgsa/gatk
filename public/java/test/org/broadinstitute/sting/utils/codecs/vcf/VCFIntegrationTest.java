package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VCFIntegrationTest extends WalkerTest {

    @Test(enabled = true)
    public void testReadingAndWritingWitHNoChanges() {

        String md5ofInputVCF = "babf02baabcfa7f72a2c6f7da5fdc996";
        String testVCF = testDir + "vcf4.1.example.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T VariantAnnotator --variant " + testVCF + " -L " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList(md5ofInputVCF));
        List<File> result = executeTest("Test Variant Annotator with no changes", spec1).getFirst();

        String test2 = baseCommand + "-T VariantsToVCF --variant " + result.get(0).getAbsolutePath();
        WalkerTestSpec spec2 = new WalkerTestSpec(test2, 1, Arrays.asList(md5ofInputVCF));
        executeTest("Test Variants To VCF from new output", spec2);
    }

    @Test
    // See https://getsatisfaction.com/gsa/topics/support_vcf_4_1_structural_variation_breakend_alleles?utm_content=topic_link&utm_medium=email&utm_source=new_topic
    public void testReadingAndWritingBreakpointAlleles() {
        String testVCF = testDir + "breakpoint-example.vcf";
        //String testVCF = validationDataLocation + "multiallelic.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("355b029487c3b4c499140d71310ca37e"));
        executeTest("Test reading and writing breakpoint VCF", spec1);
    }

    @Test
    public void testReadingAndWritingSamtools() {
        String testVCF = testDir + "samtools.vcf";

        String baseCommand = "-R " + b37KGReference + " --no_cmdline_in_header -o %s ";

        String test1 = baseCommand + "-T SelectVariants -V " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList("0807ff11bebba81b87a273ad6bee01a8"));
        executeTest("Test reading and writing samtools vcf", spec1);
    }
}
