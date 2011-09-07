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
    public void testVariantsToVCFUsingDbsnpInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("d64942fed2a5b7b407f9537dd2b4832e");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:OldDbsnp " + GATKDataLocation + "Comparisons/Validated/dbSNP/dbsnp_129_b36.rod" +
                        " -T VariantsToVCF" +
                        " -L 1:1-30,000,000" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingDbsnpInput", spec).getFirst();
    }

    @Test
    public void testVariantsToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("4accae035d271b35ee2ec58f403c68c6");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput - calls", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("2413f036ec4100b8d5db179946159a82");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,100,000-10,200,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput - genotypes", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingHapMapInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("f343085305e80c7a2493422e4eaad983");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:RawHapMap " + validationDataLocation + "rawHapMap.yri.chr1.txt" +
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
        md5.add("86f02e2e764ba35854cff2aa05a1fdd8");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:VCF " + validationDataLocation + "complexExample.vcf4" +
                        " -T VariantsToVCF" +
                        " -o %s" +
                        " -NO_HEADER",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingVCFInput", spec).getFirst();
    }
}
