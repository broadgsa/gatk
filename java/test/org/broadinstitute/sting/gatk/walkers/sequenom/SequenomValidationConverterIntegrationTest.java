package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class SequenomValidationConverterIntegrationTest extends WalkerTest {
    @Test
    public void testSNPs() {
        String testPedFile = validationDataLocation + "Sequenom_Test_File.txt";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B sequenom,Plink,"+testPedFile+" -o %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("2dab4630f40b76c0762de83fcbb60d09"));
        executeTest("Test SNPs", spec);
    }

    @Test
    public void testIndels() {
        String testPedFile = validationDataLocation + "pilot2_indel_validation.renamed.ped";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B sequenom,Plink,"+testPedFile+" -o %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("fad2dd71550dec064d458c4aa83e4de9"));
        executeTest("Test Indels", spec);
    }
}
