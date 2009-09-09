package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IntervalsTest extends WalkerTest {
    @Test
    public void testIntervals() {

        WalkerTest.WalkerTestSpec spec1 = new WalkerTest.WalkerTestSpec(
                "-T IndelIntervals -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -o %s",
                 1,
                 Arrays.asList("a4a795755b18f4ecdbc50975612bd819"));
        executeTest("testIndelIntervals", spec1);

        WalkerTest.WalkerTestSpec spec2 = new WalkerTest.WalkerTestSpec(
                "-T MismatchIntervals -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -o %s",
                 1,
                 Arrays.asList("31e8b5d4c42f2c63c08b8f6b8e10ac99"));
        executeTest("testMismatchIntervals", spec2);

        WalkerTest.WalkerTestSpec spec3 = new WalkerTest.WalkerTestSpec(
                "-T IntervalMerger -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -L 1:10,000,000-11,000,000 -intervals /humgen/gsa-scr1/GATK_Data/Validation_Data/indelIntervals.test -intervals /humgen/gsa-scr1/GATK_Data/Validation_Data/mismatchIntervals.test -o %s",
                 1,
                 Arrays.asList("bf1f23667ef0065bbcb9754f50c2d664"));
        executeTest("testMergeIntervals", spec3);

    }
}
