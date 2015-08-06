/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.coverage;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
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

    @Test
    public void testMapQ0Only() {
        String[] intervals = {"1:10,000,000-10,002,000","1:10,003,000-10,004,000"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam"};

        String cmd = buildRootCmd(b36KGReference,new ArrayList<String>(Arrays.asList(bams)), new ArrayList<String>(Arrays.asList(intervals))) + " --maxMappingQuality 0";

        WalkerTestSpec spec = new WalkerTestSpec(cmd,0,new ArrayList<String>());

        // our base file
        final File baseOutputFile = createTempFile("depthofcoveragemapq0",".tmp");

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
        final String[] intervals = {"1:1105290-1105295"};
        final String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/pilot3.CEU+TSI.5loci.bam"};
        final String cmd = buildRootCmd(b36KGReference, new ArrayList<String>(Arrays.asList(bams)), new ArrayList<String>(Arrays.asList(intervals)));

        final WalkerTestSpec spec = new WalkerTestSpec(cmd,0,new ArrayList<String>());

        final File baseOutputFile = createTempFile("testManySamples",".tmp");
        
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("c9561b52344536d2b06ab97b0bb1a234",baseOutputFile);

        execute("testLotsOfSamples",spec);
    }

}
