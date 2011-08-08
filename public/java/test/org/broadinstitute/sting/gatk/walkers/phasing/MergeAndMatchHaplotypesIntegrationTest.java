package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class MergeAndMatchHaplotypesIntegrationTest extends WalkerTest {
    private static String mergeAndMatchHaplotypesTestDataRoot = validationDataLocation + "/MergeAndMatchHaplotypes";
    private static String fundamentalTestPBTVCF = mergeAndMatchHaplotypesTestDataRoot + "/" + "FundamentalsTest.pbt.vcf";
    private static String fundamentalTestRBPVCF = mergeAndMatchHaplotypesTestDataRoot + "/" + "FundamentalsTest.pbt.rbp.vcf";

    @Test
    public void testBasicFunctionality() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T MergeAndMatchHaplotypes",
                        "-R " + b37KGReference,
                        "-B:pbt " + fundamentalTestPBTVCF,
                        "-B:rbp " + fundamentalTestRBPVCF,
                        "-o %s"
                ),
                1,
                Arrays.asList("")
        );
        executeTest("testBasicMergeAndMatchHaplotypesFunctionality", spec);
    }
}
