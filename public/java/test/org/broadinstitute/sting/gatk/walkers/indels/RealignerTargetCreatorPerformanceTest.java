package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class RealignerTargetCreatorPerformanceTest extends WalkerTest {
    @Test
    public void testRealignerTargetCreator() {

        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T RealignerTargetCreator" +
                        " --known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-50,000,000" +
                        " -o /dev/null",
                 0,
                new ArrayList<String>(0));
        executeTest("testRealignerTargetCreatorWholeGenome", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T RealignerTargetCreator" +
                        " --known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L " + evaluationDataLocation + "whole_exome_agilent_designed_120.targets.chr1.interval_list" +
                        " -o /dev/null",
                 0,
                new ArrayList<String>(0));
        executeTest("testRealignerTargetCreatorWholeExome", spec2);
    }
}
