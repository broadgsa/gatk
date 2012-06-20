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
                Arrays.asList("cd112ec37a9e28d366aff29a85fdcaa0","313cc749c7ee97713e4551de39e01ac5")
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
                Arrays.asList("27ccd6feb51de7e7dcdf35f4697fa4eb","dd90dad9fd11e1b16e6660c3ca0553e7")
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
                Arrays.asList("719d681bb0a52a40bc854bba107c5c94","b35a86d2cad17f0db7b5e84ddc0e5545")
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
                Arrays.asList("7f4a277aee2c7398fcfa84d6c98d5fb3","c53b5fd377bef48e9c6035a94db398db")
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
                Arrays.asList("44e09d2f9e4d8a9488226d03a97fe999","6f596470740e1a57679bbb38c0126364")
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
                Arrays.asList("b35a86d2cad17f0db7b5e84ddc0e5545")
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
                Arrays.asList("60ced3d078792a150a03640b62926857","6d550784382aa910f78b533d889c91c0")
        );
        executeTest("testFatherAlleleFirst", spec);
    }

}
