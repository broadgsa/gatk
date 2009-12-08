package org.broadinstitute.sting.utils.genotype.vcf;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class VCFIntegrationTest extends WalkerTest {

    @Test
    public void test1() {
        // Read in and then emit each record
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintRODs -R /broad/1KG/reference/human_b36_both.fasta -L 1:10,000,000-10,050,000 -o %s -B vcf,VCF,/humgen/gsa-scr1/GATK_Data/Validation_Data/complexExample.vcf", 1,
                Arrays.asList("26ad7a663d0f247ac26ce5490edd7ec0"));
        executeTest("test vcf", spec);
    }
}