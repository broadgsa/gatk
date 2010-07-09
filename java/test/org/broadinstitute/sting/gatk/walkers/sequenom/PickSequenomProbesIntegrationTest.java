package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class PickSequenomProbesIntegrationTest extends WalkerTest {
    @Test
    public void testProbes() {
        String testVCF = validationDataLocation + "complexExample.vcf";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T PickSequenomProbes -L 1:10,000,000-11,000,000 -B input,VCF,"+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("0f356354a4a78ff62b2848431ec11262"));
        executeTest("Test probes", spec);
    }

    @Test
    public void testProbesUsingDbSNPMask() {
        String testVCF = validationDataLocation + "pickSeqIntegrationTest.vcf";
        String testArgs = "-snp_mask " + GATKDataLocation + "/dbsnp_130_b36.rod -R "
                + oneKGLocation + "reference/human_b36_both.fasta -omitWindow -nameConvention "
                + "-project_id 1kgp3_s4_lf -T PickSequenomProbes -L " + validationDataLocation +
                "pickSeqIntegrationTest.interval_list -B input,VCF4,"+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("8b5b715b9918a0b70f4868614f197b72"));
        executeTest("Test probes", spec);
    }

    // 03c8cef968ae2d0ef5f51ac82b24f891

    @Test
    public void testProbesUsingDbSNPMaskWithNMW1() {
    String testVCF = validationDataLocation + "pickSeqIntegrationTest.vcf";
        String testArgs = "-snp_mask " + GATKDataLocation + "/dbsnp_130_b36.rod -R "
                + oneKGLocation + "reference/human_b36_both.fasta -omitWindow -nameConvention "
                + "-nmw 1 -project_id 1kgp3_s4_lf -T PickSequenomProbes -L " + validationDataLocation +
	    "pickSeqIntegrationTest.interval_list -B input,VCF4,"+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
						 Arrays.asList("03c8cef968ae2d0ef5f51ac82b24f891"));
        executeTest("Test probes", spec);
    }
}
