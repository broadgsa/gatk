package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.List;
import java.util.ArrayList;


/**
 * @author aaron
 *         <p/>
 *         Class VariantsToVCFIntegrationTest
 *         <p/>
 *         test(s) for the VariantsToVCF walker.
 */
public class VariantsToVCFIntegrationTest extends WalkerTest {


    @Test
    public void testVariantsToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("815b82fff92aab41c209eedce2d7e7d9");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:variant,GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput #1", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("22336ee9c12aa222ce29c3c5babca7d0");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:variant,GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput #2", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingHapMapInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("9bedaa7670b86a07be5191898c3727cf");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:variant,HapMap " + validationDataLocation + "rawHapMap.yri.chr1.txt" +
                        " -T VariantsToVCF" +
                        " -L 1:1-1,000,000" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingHapMapInput", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingVCFInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("cc215edec9ca28e5c79ab1b67506f9f7");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -B:variant,VCF " + validationDataLocation + "complexExample.vcf4" +
                        " -T VariantsToVCF" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingVCFInput", spec).getFirst();
    }
}
