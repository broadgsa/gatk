package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

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
        md5.add("519593d09da03e6503a863dce439151b");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,GeliText," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput #1", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingGeliInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("4541686d38eced70b8fb6647551d2329");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,GeliText," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput #2", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingHapMapInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("28728ad3a6af20a1e1aaaf185ffbff2b");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,HapMap," + validationDataLocation + "rawHapMap.yri.chr1.txt" +
                        " -T VariantsToVCF" +
                        " -L 1:1-1,000,000" +
                        " -o %s",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingHapMapInput", spec).getFirst();
    }

    @Test
    public void testGenotypesToVCFUsingVCFInput() {
        List<String> md5 = new ArrayList<String>();
        md5.add("919eb499bfcc980a14825a0265e575e3");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,VCF," + validationDataLocation + "complexExample.vcf4" +
                        " -T VariantsToVCF" +
                        " -o %s",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingVCFInput", spec).getFirst();
    }
}
