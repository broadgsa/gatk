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
                        "--no_cmdline_in_header",
                        "-R " + b37KGReference,
                        "--variant " + TNTest,
                        "-ped "+ goodFamilyFile,
                        "-L 1:10109-10315",
                        "-mvf %s",
                        "-o %s"
                ),
                2,
                Arrays.asList("d54a142d68dca54e478c13f9a0e4c95c","1a37fcc93a73429f9065b942ab771233")
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
                Arrays.asList("883ea7fd2b200c4b7fa95a4f7aa15931","7b1f5309c3d4f4aa7e9061f288dceb68")
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
                Arrays.asList("e812d62a3449b74b6948ee7deb8a0790","d00922496759e84c66a4b5e222e36997")
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
                Arrays.asList("e3c572f933a40e1878a2cfa52049517a","60e4f0be344fb944ab3378f9ab27da64")
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
                Arrays.asList("b42af3b73a2cb38cfc92f8047dd686b3","a69c3f9c005e852b44c29ab25e87ba0d")
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
                Arrays.asList("d00922496759e84c66a4b5e222e36997")
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
                Arrays.asList("c158a3816357597543ef85c4478c41e8","4f8daca19c8f31bd87850c124f91e330")
        );
        executeTest("testFatherAlleleFirst", spec);
    }

}
