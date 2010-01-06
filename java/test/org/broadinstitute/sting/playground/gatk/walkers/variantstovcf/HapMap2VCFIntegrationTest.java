package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.List;
import java.util.ArrayList;
import java.io.File;


/**
 * @author aaron
 *         <p/>
 *         Class VariantsToVCFIntegrationTest
 *         <p/>
 *         test(s) for the VariantsToVCF walker.
 */
public class HapMap2VCFIntegrationTest extends WalkerTest {


    @Test
    public void testHapMap2VCF() {
        List<String> md5 = new ArrayList<String>();
        md5.add("4d36df142bbd3d446baec6213771800a");

        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " --rodBind HapMapChip,HapMapGenotype," + validationDataLocation + "hapmap_phase_ii+iii_genotypes_chrX_YRI_r27_nr.hapmap" +
                        " -T HapMap2VCF" +
                        " -L chrX:1-1,000,000" +
                        " --vcfOutput %s",
                1, // just one output file
                md5);
        List<File> result = executeTest("testHapMap2VCFUsingGeliInput", spec).getFirst();
    }

}