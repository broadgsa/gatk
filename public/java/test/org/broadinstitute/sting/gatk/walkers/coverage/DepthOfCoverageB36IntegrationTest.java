package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

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
        String[] intervals = {"1:10,000,000-10,002,000","1:10,003,000-10,004,000"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam"};

        String cmd = buildRootCmd(b36KGReference,new ArrayList<String>(Arrays.asList(bams)), new ArrayList<String>(Arrays.asList(intervals))) + " --maxMappingQuality 0";

        WalkerTestSpec spec = new WalkerTestSpec(cmd,0,new ArrayList<String>());

        // our base file
        File baseOutputFile = this.createTempFile("depthofcoveragemapq0",".tmp");

        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("f39af6ad99520fd4fb27b409ab0344a0",baseOutputFile);
        spec.addAuxFile("6b15f5330414b6d4e2f6caea42139fa1", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("cc6640d82077991dde8a2b523935cdff", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("0fb627234599c258a3fee1b2703e164a", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("cb73a0fa0cee50f1fb8f249315d38128", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("347b47ef73fbd4e277704ddbd7834f69", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("4ec920335d4b9573f695c39d62748089", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));

        execute("testMapQ0Only",spec);
    }

    @Test
    public void testLotsOfSamples() {
        String[] intervals = {"1:1105290-1105295"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/pilot3.CEU+TSI.5loci.bam"};
        String cmd = buildRootCmd(b36KGReference, new ArrayList<String>(Arrays.asList(bams)), new ArrayList<String>(Arrays.asList(intervals)));

        WalkerTestSpec spec = new WalkerTestSpec(cmd,0,new ArrayList<String>());

        File baseOutputFile = this.createTempFile("testManySamples",".tmp");
        
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("c9561b52344536d2b06ab97b0bb1a234",baseOutputFile);

        execute("testLotsOfSamples",spec);
    }

}
