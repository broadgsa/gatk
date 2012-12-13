package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class PhaseByTransmissionIntegrationTest extends WalkerTest {
    private static String phaseByTransmissionTestDataRoot = privateTestDir + "PhaseByTransmission/";
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
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + TNTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("af979bcb353edda8dee2127605c71daf","1ea9994f937012e8de599ec7bcd62a0e")
        );
        executeTest("testTrueNegativeMV", spec);
    }

    @Test
    public void testTruePositiveMV() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + TPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("1dc36ff8d1d5f5d2c1c1bf21517263bf","547fdfef393f3045a96d245ef6af8acb")
        );
        executeTest("testTruePositiveMV", spec);
    }

    @Test
    public void testFalsePositiveMV() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("ae60f2db6102ca1f4e93cd18d0634d7a","9529e2bf214d72e792d93fbea22a3b91")
        );
        executeTest("testFalsePositiveMV", spec);
    }

    @Test
    public void testSpecialCases() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + SpecialTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("590ee56e745984296f73e4277277eac7","8c157d79dd00063d2932f0d2b96f53d8")
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
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-prior 1e-4",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("78158d738917b8f0b7a736a1739b2cc5","343e418850ae4a687ebef2acd55fcb07")
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
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + FPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-o %s"
                ),
                1,
                Arrays.asList("9529e2bf214d72e792d93fbea22a3b91")
        );
        executeTest("testMVFileOption", spec);
    }

    //Test when running with the fatherAlleleFirst option
    @Test
    public void testFatherAlleleFirst() {
        WalkerTestSpec spec = new WalkerTestSpec(
                buildCommandLine(
                        "-T PhaseByTransmission",
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + TPTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s",
                        "-fatherAlleleFirst"
                ),
                2,
                Arrays.asList("dc6afb769b55e6038677fa590b2b2e89","52ffa82428e63ade22ea37b72ae58492")
        );
        executeTest("testFatherAlleleFirst", spec);
    }

}
