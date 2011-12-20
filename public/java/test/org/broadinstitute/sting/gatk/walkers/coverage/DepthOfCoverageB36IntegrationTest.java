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
        spec.addAuxFile("5b6c16a1c667c844882e9dce71454fc4",baseOutputFile);
        spec.addAuxFile("fc161ec1b61dc67bc6a5ce36cb2d02c9", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("89321bbfb76a4e1edc0905d50503ba1f", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("0fb627234599c258a3fee1b2703e164a", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("4dd16b659065e331ed4bd3ab0dae6c1b", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("2be0c18b501f4a3d8c5e5f99738b4713", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("5a26ef61f586f58310812580ce842462", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));


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
        spec.addAuxFile("d73fa1fc492f7dcc1d75056f8c12c92a",baseOutputFile);

        execute("testLotsOfSamples",spec);
    }

}
