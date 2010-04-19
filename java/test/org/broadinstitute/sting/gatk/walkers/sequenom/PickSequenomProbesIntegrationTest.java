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
                Arrays.asList("6b5409cc78960f1be855536ed89ea9dd"));
        executeTest("Test probes", spec);
    }
}