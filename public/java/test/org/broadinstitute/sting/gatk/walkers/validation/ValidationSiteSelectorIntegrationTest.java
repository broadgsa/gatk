package org.broadinstitute.sting.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: 3/26/12
 * Time: 3:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class ValidationSiteSelectorIntegrationTest extends WalkerTest {
    public static String baseTestString(String args) {
        return "-T ValidationSiteSelector -R " + b36KGReference + " -L 1 -o %s --no_cmdline_in_header -numSites 100 " + args;
    }

    private static String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
    private static String samplesFile = validationDataLocation + "SelectVariants.samples.txt";
    private static String samplePrefix = " -sf " + samplesFile;
    private static String freqUnif = " --frequencySelectionMode UNIFORM ";
    private static String freqAF = " --frequencySelectionMode KEEP_AF_SPECTRUM ";
    private static String sampleNone = " -sampleMode NONE ";
    private static String sampleGT = samplePrefix+" -sampleMode POLY_BASED_ON_GT ";
    private static String sampleGL = samplePrefix+" -sampleMode POLY_BASED_ON_GL -samplePNonref 0.95";


    @Test(enabled=true)
    public void testNoSampleSelectionFreqUniform() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleNone + freqUnif + "--variant " + testfile),
                1,
                Arrays.asList("b8a988757ac1f206d123140da5a3e778")
        );

        executeTest("testNoSampleSelectionFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testNoSampleSelectionFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleNone + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("542d5d5ff8c64da7b077bab4b950a9a3")
        );

        executeTest("testNoSampleSelectionFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqUniform() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqUnif + "--variant " + testfile),
                1,
                Arrays.asList("7385b17eed7f4ff0f6e82e60c3334ce7")
        );

        executeTest("testPolyGTFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("0ee4a565a0d4f6b6942abd72a373becd")
        );

        executeTest("testPolyGTFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGLFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGL + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("3bf094e1aef563daf7c936032259d490")
        );

        executeTest("testPolyGLFreqAF--" + testfile, spec);
    }

}
