package org.broadinstitute.sting.oneoffprojects.walkers.coverage;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 25, 2010
 */
public class CoverageStatisticsIntegrationTest extends WalkerTest {

    private boolean RUN_TESTS = true;
    private String root = "-T CoverageStatistics ";

    private String buildRootCmd(String ref, String bam, String interval) {
        return root + "-R "+ref+" -I "+bam+" -L "+interval;
    }

    private void execute(String name, WalkerTestSpec spec) {
        if ( RUN_TESTS ) {
            executeTest(name,spec);
        }
    }

    @Test
    public void testBaseOutputNoFiltering() {
        // our base file
        File baseOutputFile = this.createTempFile("outputtemp",".tmp");
        this.setOutputFileLocation(baseOutputFile);

        String bam_file = "/humgen/gsa-hphome1/chartl/projects/depthOfCoverage/testFiles/bams/Ciliopathies_1_88534_3_samples.bam";
        String interval_list = "chr1:855534";
        String reference = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
        String cmd = buildRootCmd(reference,bam_file,interval_list) + " -mmq 0 -mbq 0 -omitSampleSummary -omitLocus";
        String expected = "d41d8cd98f00b204e9800998ecf8427e";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,1, Arrays.asList(expected));

        // now add the expected files that get generated
        spec.addAuxFile("344936e0bb4613544ff137cd1a002d6c",new File(baseOutputFile.getAbsolutePath() + ".interval_statistics"));

        execute("testBaseOutputNoFiltering",spec);
    }
}
