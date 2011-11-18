package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PhaseByTransmissionIntegrationTest extends WalkerTest {
    private static String phaseByTransmissionTestDataRoot = validationDataLocation + "PhaseByTransmission/";
    private static String goodFamilyFile =  phaseByTransmissionTestDataRoot + "PhaseByTransmission.IntegrationTest.goodFamilies.ped";
    private static String TNTest = phaseByTransmissionTestDataRoot + "PhaseByTransmission.IntegrationTest.TN.vcf";
    private static String TPTest = phaseByTransmissionTestDataRoot + "PhaseByTransmission.IntegrationTest.TP.vcf";
    private static String FPTest = phaseByTransmissionTestDataRoot + "PhaseByTransmission.IntegrationTest.FP.vcf";
    private static String SpecialTest = phaseByTransmissionTestDataRoot + "PhaseByTransmission.IntegrationTest.Special.vcf";

    //Tests using PbT on all genotypes with default parameters
    //And all reporting options
    @Test
    public void testTrueNegativeMV() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + TNTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("16fefda693156eadf1481fd9de23facb","9418a7a6405b78179ca13a67b8bfcc14")
        );
        executeTest("testTrueNegativeMV", spec);
    }

    @Test
    public void testTruePositiveMV() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + TPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("14cf1d21a54d8b9fb506df178b634c56","efc66ae3d036715b721f9bd35b65d556")
        );
        executeTest("testTruePositiveMV", spec);
    }

    @Test
    public void testFalsePositiveMV() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("f9b0fae9fe1e0f09b883a292b0e70a12","398724bc1e65314cc5ee92706e05a3ee")
        );
        executeTest("testFalsePositiveMV", spec);
    }

    @Test
    public void testSpecialCases() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + SpecialTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("b8d1aa3789ce77b45430c62d13ee3006","a1a333e08fafb288cda0e7711909e1c3")
        );
        executeTest("testSpecialCases", spec);
    }

    //Test using a different prior
    //Here the FP file is used but as the prior is lowered, 3 turn to TP
    @Test
    public void testPriorOption() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-prior 1e-4",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("7201ce7cc47db5840ac6b647709f7c33","c11b5e7cd7459d90d0160f917eff3b1e")
        );
        executeTest("testPriorOption", spec);
    }

    //Test when running without MV reporting option
    //This is the exact same test file as FP but should not generate a .mvf file
    @Test
    public void testMVFileOption() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "-NO_HEADER",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-o %s"
                ),
                1,
                Arrays.asList("398724bc1e65314cc5ee92706e05a3ee")
        );
        executeTest("testMVFileOption", spec);
    }

}
