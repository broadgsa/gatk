package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VCFIntegrationTest extends WalkerTest {

    @Test
    public void test1() {
        // Read in and then emit each record
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintRODs -R " + oneKGLocation + "reference/human_b36_both.fasta -L 1:10,000,000-10,050,000 -o %s -B vcf,VCF," + validationDataLocation + "complexExample.vcf", 1,
                Arrays.asList("68b123acca4975553297fcd776c70464"));
        executeTest("test vcf", spec);
    }
}