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
        return "-T ValidationSiteSelector -R " + b36KGReference + " -L 1 -o %s -NO_HEADER -numSites 100 " + args;
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
                Arrays.asList("d49baeb8000a426c172ce1d81eb37963")
        );

        executeTest("testNoSampleSelectionFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testNoSampleSelectionFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleNone + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("0fb0d015d462c34514fc7e96beea5f56")
        );

        executeTest("testNoSampleSelectionFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqUniform() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqUnif + "--variant " + testfile),
                1,
                Arrays.asList("0672854299d42ea8af906976a3849ae6")
        );

        executeTest("testPolyGTFreqUniform--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGTFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGT + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("5bdffda1a063d0bddd6b236854ec627d")
        );

        executeTest("testPolyGTFreqAF--" + testfile, spec);
    }

    @Test(enabled=true)
    public void testPolyGLFreqAF() {

        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString(sampleGL + freqAF + "--variant " + testfile),
                1,
                Arrays.asList("35ef16aa41303606a4b94f7b88bd9aa8")
        );

        executeTest("testPolyGLFreqAF--" + testfile, spec);
    }

}
