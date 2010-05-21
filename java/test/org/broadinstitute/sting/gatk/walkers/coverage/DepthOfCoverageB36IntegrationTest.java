package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date May 20, 2010
 */
public class DepthOfCoverageB36IntegrationTest extends WalkerTest {

    private boolean RUN_TESTS = true;
    private String root = "-T DepthOfCoverage ";
    private String hg18 = "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    private String b36 = "/broad/1KG/reference/human_b36_both.fasta";

    private String buildRootCmd(String ref, List<String> bams, List<String> intervals) {
        StringBuilder bamBuilder = new StringBuilder();
        do {
            bamBuilder.append(" -I ");
            bamBuilder.append(bams.remove(0));
        } while ( bams.size() > 0 );

        StringBuilder intervalBuilder = new StringBuilder();
        do {
            intervalBuilder.append(" -L ");
            intervalBuilder.append(intervals.remove(0));
        } while ( intervals.size() > 0 );


        return root + "-R "+ref+bamBuilder.toString()+intervalBuilder.toString();
    }

    private void execute(String name, WalkerTest.WalkerTestSpec spec) {
        if ( RUN_TESTS ) {
            executeTest(name,spec);
        }
    }

    public File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }

    @Test
    public void testMapQ0Only() {
        // our base file
        File baseOutputFile = this.createTempFile("depthofcoveragemapq0",".tmp");
        this.setOutputFileLocation(baseOutputFile);

        String[] intervals = {"1:10,000,000-10,002,000","1:10,003,000-10,004,000"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam"};

        String cmd = buildRootCmd(b36,new ArrayList<String>(Arrays.asList(bams)), new ArrayList<String>(Arrays.asList(intervals))) + " --maxMappingQuality 0";

        WalkerTestSpec spec = new WalkerTestSpec(cmd,0,new ArrayList<String>());

        spec.addAuxFile("f39af6ad99520fd4fb27b409ab0344a0",baseOutputFile);
        spec.addAuxFile("6b15f5330414b6d4e2f6caea42139fa1", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("cc6640d82077991dde8a2b523935cdff", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("0fb627234599c258a3fee1b2703e164a", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("6f6085293e8cca402ef4fe648cf06ea8", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("347b47ef73fbd4e277704ddbd7834f69", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("ff0c79d161259402d54d572151582697", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));

        execute("testMapQ0Only",spec);
    }

}
