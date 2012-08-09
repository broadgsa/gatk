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
                Arrays.asList("f4b0b5471e03306ee2fad27d88b217b6","f8721f4f5d3bae2848ae15c3f120709b")
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
                Arrays.asList("dbc64776dcc9e01a468b61e4e0db8277","547fdfef393f3045a96d245ef6af8acb")
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
                Arrays.asList("37793e78861bb0bc070884da67dc10e6","9529e2bf214d72e792d93fbea22a3b91")
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
                Arrays.asList("e4da7639bb542d6440975da12b94973f","8c157d79dd00063d2932f0d2b96f53d8")
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
                Arrays.asList("ab92b714471a000285577d540e1fdc2e","343e418850ae4a687ebef2acd55fcb07")
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
                Arrays.asList("4b937c1b4e96602a7479b07b59254d06","52ffa82428e63ade22ea37b72ae58492")
        );
        executeTest("testFatherAlleleFirst", spec);
    }

}
