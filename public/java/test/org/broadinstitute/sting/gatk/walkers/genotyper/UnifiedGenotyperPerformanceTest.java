package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class UnifiedGenotyperPerformanceTest extends WalkerTest {

    @Test
    public void testUnifiedGenotyperWholeGenome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T UnifiedGenotyper" +
                        " -glm BOTH" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " --dbsnp:VCF " + hg18dbSNP132 +
                        " -o /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testUnifiedGenotyperWholeGenome", spec);
    }

    @Test
    public void testUnifiedGenotyperWholeExome() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T UnifiedGenotyper" +
                        " -glm BOTH" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " --dbsnp:vcf " + hg18dbSNP132 +
                        " -o /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testUnifiedGenotyperWholeExome", spec);
    }

    @Test
    public void testUnifiedGenotyperWGParallel() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T UnifiedGenotyper" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -glm BOTH" +
                        " -L chr1:1-50,000,000" +
                        " -nt 10" +
                        " --dbsnp:vcf " + hg18dbSNP132 +
                        " -o /dev/null",
                0,
                new ArrayList<String>(0));
        executeTest("testUnifiedGenotyperWGParallel", spec);
    }

}
