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
                 Arrays.asList("49431d567524d6fd32a569504b25f212"));
        executeTest("testBamToFasta", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T BamToFastq -reverse -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,100-10,000,500;1:10,100,000-10,101,000;1:10,900,000-10,900,001 -o %s",
                 1,
                 Arrays.asList("f3a4a39d36270136c12bb1315fdb7dff"));
        executeTest("testBamToFastaReverse", spec2);
    }
}
