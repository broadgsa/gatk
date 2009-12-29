package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class BamToFastqIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T BamToFastq -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,100-10,000,500;1:10,100,000-10,101,000;1:10,900,000-10,900,001 -o %s",
                 1,
                 Arrays.asList("2ec2ff5ac405099bf48186d2c7e5fabd"));
        executeTest("testBamToFasta", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T BamToFastq -reverse -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,100-10,000,500;1:10,100,000-10,101,000;1:10,900,000-10,900,001 -o %s",
                 1,
                 Arrays.asList("5a79484da43925b0e0461be28fdad07c"));
        executeTest("testBamToFastaReverse", spec2);
    }
}
