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
                Arrays.asList("6a9e990a9252840904b5144213915b32")
        );

        executeTest("testNoSampleSelectionFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testNoSampleSelectionFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleNone + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("eaa2385086cddff68cf4fdb81cbdbbb9")
        );

        executeTest("testNoSampleSelectionFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqUniform() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqUnif + "--variant " + testfile),
                1,
                Arrays.asList("24077656f590d6905546f7e019c8dccb")
        );

        executeTest("testPolyGTFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("3c1180fd9b5e80e540b39c5a95fbe722")
        );

        executeTest("testPolyGTFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGLFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGL + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("ad30c028864348204ebe80b9c8c503e8")
        );

        executeTest("testPolyGLFreqAF--" + testfile, spec);
    }

}
