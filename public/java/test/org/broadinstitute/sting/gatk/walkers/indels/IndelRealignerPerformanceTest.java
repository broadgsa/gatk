package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.ArrayList;

public class IndelRealignerPerformanceTest extends WalkerTest {
    @Test
    public void testHighCoverage() {
        WalkerTestSpec spec = new WalkerTestSpec(

                "-R " + b36KGReference +
                        " -T IndelRealigner" +
                        " -I " + validationDataLocation + "indelRealignerTest.pilot1.veryHighCoverage.bam" +
                        " -L 20:49,500-55,500" +
                        " -o /dev/null" +
                        " -targetIntervals " + validationDataLocation + "indelRealignerTest.pilot1.ceu.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerHighCoverage", spec);
    }

    @Test
    public void testRealigner() {
        WalkerTestSpec spec1 = new WalkerTestSpec(

                "-R " + hg18Reference +
                        " -T IndelRealigner" +
                        " -LOD 5" +
                        " -maxConsensuses 100" +
                        " -greedy 100" +
                        " -known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -o /dev/null" +
                        " -I " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.bam" +
                        " -L chr1:1-5,650,000" +
                        " -targetIntervals " + evaluationDataLocation + "NA12878.GAII.chr1.50MB.realigner.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerWholeGenome", spec1);

        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T IndelRealigner" +
                        " -LOD 5" +
                        " -maxConsensuses 100" +
                        " -greedy 100" +
                        " -known " + GATKDataLocation + "dbsnp_132.hg18.vcf" +
                        " -o /dev/null" +
                        " -I " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.bam" +
                        " -L chr1:1-150,000,000" +
                        " -targetIntervals " + evaluationDataLocation + "NA12878.ESP.WEx.chr1.realigner.intervals",
                 0,
                new ArrayList<String>(0));
        executeTest("testIndelRealignerWholeExome", spec2);
    }
}
