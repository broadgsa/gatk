package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class VCFIntegrationTest extends WalkerTest {

    @Test
    public void testReadingAndWritingWitHNoChanges() {

        String md5ofInputVCF = "a990ba187a69ca44cb9bc2bb44d00447";
        String testVCF = validationDataLocation + "vcf4.1.example.vcf";

        String baseCommand = "-R " + b37KGReference + " -NO_HEADER -o %s ";

        String test1 = baseCommand + "-T VariantAnnotator --variant " + testVCF + " -L " + testVCF;
        WalkerTestSpec spec1 = new WalkerTestSpec(test1, 1, Arrays.asList(md5ofInputVCF));
        List<File> result = executeTest("Test Variant Annotator with no changes", spec1).getFirst();

        String test2 = baseCommand + "-T VariantsToVCF --variant " + result.get(0).getAbsolutePath();
        WalkerTestSpec spec2 = new WalkerTestSpec(test2, 1, Arrays.asList(md5ofInputVCF));
        executeTest("Test Variants To VCF from new output", spec2);
    }
}
