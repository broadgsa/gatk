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
                Arrays.asList("2e273d400b4b69e39c34e465b200b192"));
        executeTest("Test SNPs", spec);
    }

    @Test
    public void testIndels() {
        String testPedFile = validationDataLocation + "pilot2_indel_validation.renamed.ped";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B sequenom,Plink,"+testPedFile+" -o %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("e15a63fc49ec25ebcae60a28a5f3f830"));
        executeTest("Test Indels", spec);
    }
}
