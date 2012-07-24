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
        md5.add("72e6ce7aff7dec7ca9e7580be7ddd435");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:OldDbsnp " + GATKDataLocation + "Comparisons/Validated/dbSNP/dbsnp_129_b36.rod" +
                        " -T VariantsToVCF" +
                        " -L 1:1-30,000,000" +
                        " -o %s" +
                        " --no_cmdline_in_header",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingDbsnpInput", spec).getFirst();
    }

    @Test
    public void testVariantsToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("22373883afa2221b5a4f75a50f30f26b");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " --no_cmdline_in_header",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput - calls", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("738eb66dbc400dcd1786cd9e49902e8c");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:GeliText " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,100,000-10,200,000" +
                        " -sample NA123AB" +
                        " -o %s" +
                        " --no_cmdline_in_header",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput - genotypes", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingHapMapInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("67656672acc264156f5a3e01e5cac61a");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:RawHapMap " + validationDataLocation + "rawHapMap.yri.chr1.txt" +
                        " -T VariantsToVCF" +
                        " -L 1:1-1,000,000" +
                        " -o %s" +
                        " --no_cmdline_in_header",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingHapMapInput", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingVCFInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("21084d32ce7ac5df3cee1730bfaaf71c");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " --variant:VCF " + privateTestDir + "complexExample.vcf4" +
                        " -T VariantsToVCF" +
                        " -o %s" +
                        " --no_cmdline_in_header",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingVCFInput", spec).getFirst();
    }
}
