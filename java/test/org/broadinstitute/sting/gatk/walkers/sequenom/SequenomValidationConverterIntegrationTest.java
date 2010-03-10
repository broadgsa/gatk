package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class SequenomValidationConverterIntegrationTest extends WalkerTest {
    @Test
    public void testSNPs() {
        String testPedFile = validationDataLocation + "Sequenom_Test_File.txt";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B input,Plink,"+testPedFile+" -vcf %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("d19f28fdbe3e731522a52c5329777a9f"));
        executeTest("Test SNPs", spec);
    }

    @Test
    public void testIndels() {
        String testPedFile = validationDataLocation + "pilot2_indel_validation.renamed.ped";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B input,Plink,"+testPedFile+" -vcf %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("257fcd5e345f2853813e37b88fbc707c"));
        executeTest("Test Indels", spec);
    }
}
