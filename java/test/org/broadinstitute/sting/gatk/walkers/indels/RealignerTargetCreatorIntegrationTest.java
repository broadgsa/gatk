package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class RealignerTargetCreatorIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000 -o %s --realignReadsWithBadMates",
                 1,
                 Arrays.asList("d21e83a8b0d3f63acd9ca3b0b636e515"));
        executeTest("test standard", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -D /humgen/gsa-hpprojects/GATK/data/dbsnp_129_b36.rod -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,050,000 -o %s --realignReadsWithBadMates",
                 1,
                 Arrays.asList("bfccfa50f62d10ee2fe8cfa68fb70002"));
        executeTest("test dbsnp", spec2);

        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + oneKGLocation + "reference/human_b36_both.fasta -B indels,VCF," + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf -BTI indels -o %s",
                 1,
                 Arrays.asList("1a11cfc9cc713617c82bdec503ebe02a"));
        executeTest("test rods only", spec3);
    }
}
