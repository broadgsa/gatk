package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;


public class RecalibrationWalkersPerformanceTest extends WalkerTest {

    private void testCountCovariatesWholeGenomeRunner(String moreArgs) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T CountCovariates" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -standard" +
                        " -OQ" +
                        " -knownSites " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -recalFile /dev/null" + moreArgs,
                0,
                new ArrayList<String>(0));
        executeTest("testCountCovariatesWholeGenome", spec);
    }

    private  void testCountCovariatesWholeExomeRunner(String moreArgs) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T CountCovariates" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -standard" +
                        " -OQ" +
                        " -knownSites " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -recalFile /dev/null" + moreArgs,
                0,
                new ArrayList<String>(0));
        executeTest("testCountCovariatesWholeExome", spec);
    }

    @Test
    public void testCountCovariatesWholeGenome() { testCountCovariatesWholeGenomeRunner(""); }
    @Test
    public void testCountCovariatesWholeGenomeParallel() { testCountCovariatesWholeGenomeRunner(" -nt 4"); }

    @Test
    public void testCountCovariatesWholeExome() { testCountCovariatesWholeExomeRunner(""); }
    @Test
    public void testCountCovariatesWholeExomeParallel() { testCountCovariatesWholeExomeRunner(" -nt 4"); }

    @Test
    public void testTableRecalibratorWholeGenome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T TableRecalibration" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -OQ" +
                        " -recalFile " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.recal.csv" +
                        " -o /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testTableRecalibratorWholeGenome", spec);
    }

    @Test
    public void testTableRecalibratorWholeExome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T TableRecalibration" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -OQ" +
                        " -recalFile " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.recal.csv" +
                        " -o /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testTableRecalibratorWholeExome", spec);
    }
}
