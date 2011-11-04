package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class RealignerTargetCreatorIntegrationTest extends WalkerTest {

    @Test
    public void testIntervals1() {
        String md5 = "3f0b63a393104d0c4158c7d1538153b8";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam --mismatchFraction 0.15 -L 1:10,000,000-10,050,000 -o %s",
                 1,
                 Arrays.asList(md5));
        executeTest("test standard nt=1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-nt 4 -T RealignerTargetCreator -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam --mismatchFraction 0.15 -L 1:10,000,000-10,050,000 -o %s",
                 1,
                 Arrays.asList(md5));
        executeTest("test standard nt=4", spec2);
    }

    @Test
    public void testIntervals2() {
        String md5 = "e0f745b79b679c225314a2abef4919ff";

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator --known " + b36dbSNP129 + " -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,200,000 -o %s",
                 1,
                 Arrays.asList(md5));
        executeTest("test with dbsnp nt=1", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-nt 4 -T RealignerTargetCreator --known " + b36dbSNP129 + " -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-10,200,000 -o %s",
                 1,
                 Arrays.asList(md5));
        executeTest("test with dbsnp nt=4", spec2);
    }

    @Test
    public void testKnownsOnly() {
        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                "-T RealignerTargetCreator -R " + b36KGReference + " --known " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4 -L " + validationDataLocation + "NA12878.chr1_10mb_11mb.slx.indels.vcf4 -o %s",
                 1,
                 Arrays.asList("5206cee6c01b299417bf2feeb8b3dc96"));
        executeTest("test rods only", spec3);
    }
}
