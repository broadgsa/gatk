package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class RealignerTargetCreatorIntegrationTest extends WalkerTest {

    @Test
    public void testIntervals() {

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000 -o %s",
                 1,
                 Arrays.asList("e7accfa58415d6da80383953b1a3a986"));
        executeTest("test standard", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_b36.rod -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000 -o %s",
                 1,
                 Arrays.asList("f23ba17ee0f9573dd307708175d90cd2"));
        executeTest("test dbsnp", spec2);

        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b36KGReference + " -B:indels,VCF " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4 -BTI indels -o %s",
                 1,
                 Arrays.asList("5206cee6c01b299417bf2feeb8b3dc96"));
        executeTest("test rods only", spec3);
    }
}
