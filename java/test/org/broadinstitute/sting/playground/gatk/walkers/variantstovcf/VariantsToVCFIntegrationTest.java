package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

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
        md5.add("8f32f0efed5d0233cf9292198f4f01d8");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.variants.geli.calls" +
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
        md5.add("40cc4d04d9a50043ce1322ea2650c453");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + oneKGLocation + "reference/human_b36_both.fasta" +
                        " -B variant,Variants," + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.lod5.genotypes.geli.calls" +
                        " -T VariantsToVCF" +
                        " -L 1:10,000,000-11,000,000" +
                        " -sample NA123AB" +
                        " -o %s",
                1, // just one output file
                md5);
        executeTest("testVariantsToVCFUsingGeliInput #2", spec).getFirst();
    }

}
