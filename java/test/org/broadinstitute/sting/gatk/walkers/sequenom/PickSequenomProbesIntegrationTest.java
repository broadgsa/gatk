package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PickSequenomProbesIntegrationTest extends WalkerTest {
    @Test
    public void testProbes() {
        String testVCF = validationDataLocation + "complexExample.vcf4";
        String testArgs = "-R " + b36KGReference + " -T PickSequenomProbes -L 1:10,000,000-11,000,000 -B:input,VCF "+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("6b5409cc78960f1be855536ed89ea9dd"));
        executeTest("Test probes", spec);
    }

    @Test
    public void testProbesUsingDbSNPMask() {
        String testVCF = validationDataLocation + "pickSeqIntegrationTest.vcf";
        String testArgs = "-snp_mask " + GATKDataLocation + "/dbsnp_130_b36.rod -R "
                + b36KGReference + " -omitWindow -nameConvention "
                + "-project_id 1kgp3_s4_lf -T PickSequenomProbes -L " + validationDataLocation +
                "pickSeqIntegrationTest.interval_list -B:input,VCF "+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("e50162490d0597c8f00d2d20b1463661"));
        executeTest("Test probes", spec);
    }

    @Test
    public void testProbesUsingDbSNPMaskWithNMW1() {
    String testVCF = validationDataLocation + "pickSeqIntegrationTest.vcf";
        String testArgs = "-snp_mask " + GATKDataLocation + "/dbsnp_130_b36.rod -R "
                + b36KGReference + " -omitWindow -nameConvention "
                + "-nmw 1 -project_id 1kgp3_s4_lf -T PickSequenomProbes -L " + validationDataLocation +
	    "pickSeqIntegrationTest.interval_list -B:input,VCF "+testVCF+" -o %s";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
						 Arrays.asList("e50162490d0597c8f00d2d20b1463661"));
        executeTest("Test probes", spec);
    }
}
