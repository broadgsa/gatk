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
                Arrays.asList("cb1f57e8bcaec4b599be075b6d5288a1"));
        executeTest("Test probes", spec);
    }
}
