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

        String md5 = "46d53491af1d3aa0ee1f1e13d68b732d";
        String testVCF = validationDataLocation + "pickSeqIntegrationTest.vcf";

        String testArgs = "-snp_mask " + validationDataLocation + "pickSeqIntegrationTest.bed -R "
                + b36KGReference + " -omitWindow -nameConvention "
                + "-project_id 1kgp3_s4_lf -T PickSequenomProbes -B:input,VCF "+testVCF+" -o %s";
        WalkerTestSpec spec1 = new WalkerTestSpec(testArgs, 1, Arrays.asList(md5));
        executeTest("Test probes", spec1);

        testArgs += " -nmw 1";
        WalkerTestSpec spec2 = new WalkerTestSpec(testArgs, 1, Arrays.asList(md5));
        executeTest("Test probes", spec2);
    }
}
