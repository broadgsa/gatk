package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PhaseByTransmissionIntegrationTest extends WalkerTest {
    private static String phaseByTransmissionTestDataRoot = validationDataLocation + "/PhaseByTransmission";
    private static String fundamentalTestVCF = phaseByTransmissionTestDataRoot + "/" + "FundamentalsTest.unfiltered.vcf";

    @Test
    public void testBasicFunctionalityWithoutFilters() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-R " + b37KGReference,
                        "-B:variant,VCF " + fundamentalTestVCF,
                        "-f NA12892+NA12891=NA12878",
                        "-nofilters",
                        "-o %s"
                ),
                1,
                Arrays.asList("416a483e87358cdcb0b09a496e3254c0")
        );
        executeTest("testBasicFunctionalityWithoutFilters", spec);
    }

    /*
    @Test
    public void testBasicFunctionalityWithFilters() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-R " + b37KGReference,
                        "-B:variant,VCF " + fundamentalTestVCF,
                        "-f NA12892+NA12891=NA12878",
                        "-o %s"
                ),
                1,
                Arrays.asList("8c5db343567e90e97993912c7e541d0d")
        );
        executeTest("testBasicFunctionalityWithFilters", spec);
    }
    */
}
