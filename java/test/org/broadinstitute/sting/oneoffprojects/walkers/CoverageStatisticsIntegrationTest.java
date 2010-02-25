package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 25, 2010
 */
public class CoverageStatisticsIntegrationTest extends WalkerTest {

    private boolean RUN_TESTS = false;
    private String root = "-T CoverageStatistics ";

    private String buildRootCmd(String ref, String bam, String interval) {
        return root + "-R "+ref+" -I "+bam+" -L "+interval+" -o %s";
    }

    private void execute(String name, WalkerTestSpec spec) {
        if ( RUN_TESTS ) {
            executeTest(name,spec);
        }
    }

    @Test
    public void testBaseOutputNoFiltering() {
        String bam_file = "/humgen/gsa-hphome1/chartl/projects/depthOfCoverage/testFiles/bams/Ciliopathies_1_88534_3_samples.bam";
        String interval_list = "chr1:855534";
        String reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
        String cmd = buildRootCmd(reference,bam_file,interval_list) + " -mmq 0 -mbq 0 -omitSampleSummary -omitIntervals -omitLocus";
        String expected = "2aee1dbcb69bf1e874d56cd23336afa8";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,1, Arrays.asList(expected));
        execute("testBaseOutputNoFiltering",spec);
    }
}
