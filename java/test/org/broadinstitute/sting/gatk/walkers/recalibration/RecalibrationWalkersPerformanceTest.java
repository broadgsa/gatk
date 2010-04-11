package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.ArrayList;

public class RecalibrationWalkersPerformanceTest extends WalkerTest {

    @Test
    public void testCountCovariatesWholeGenome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T CountCovariates" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -standard" +
                        " -OQ" +
                        " -recalFile /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testCountCovariatesWholeGenome", spec);
    }

    @Test
    public void testCountCovariatesWholeExome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T CountCovariates" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -standard" +
                        " -OQ" +
                        " -recalFile /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testCountCovariatesWholeExome", spec);
    }

    @Test
    public void testTableRecalibratorWholeGenome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T TableRecalibration" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -OQ" +
                        " -recalFile " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.recal.csv" +
                        " -outputBam /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testTableRecalibratorWholeGenome", spec);
    }

    @Test
    public void testTableRecalibratorWholeExome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T TableRecalibration" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -OQ" +
                        " -recalFile " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.recal.csv" +
                        " -outputBam /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testTableRecalibratorWholeExome", spec);
    }
}