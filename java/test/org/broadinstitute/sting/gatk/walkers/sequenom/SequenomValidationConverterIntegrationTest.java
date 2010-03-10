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
                Arrays.asList("02dcb683db503599efed0b76daa8fcba"));
        executeTest("Test SNPs", spec);
    }

    @Test
    public void testIndels() {
        String testPedFile = validationDataLocation + "pilot2_indel_validation.renamed.ped";
        String testArgs = "-R "+oneKGLocation+"reference/human_b36_both.fasta -T SequenomValidationConverter -B input,Plink,"+testPedFile+" -vcf %s";
        WalkerTest.WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("e72c0ab95b279a4d39cc14d40770a801"));
        executeTest("Test Indels", spec);
    }
}
