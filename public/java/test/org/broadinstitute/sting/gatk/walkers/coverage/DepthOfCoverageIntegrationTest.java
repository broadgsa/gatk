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
 * @Date Feb 25, 2010
 */
public class DepthOfCoverageIntegrationTest extends WalkerTest {

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

    private void execute(String name, WalkerTestSpec spec) {
        if ( RUN_TESTS ) {
            executeTest(name,spec);
        }
    }

    @Test
    public void testBaseOutputNoFiltering() {
        String[] intervals = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/fhs_jhs_30_targts.interval_list"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/FHS_indexed_subset.bam"};

        String cmd = buildRootCmd(hg18Reference,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 0 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv -ct 10 -ct 15 -ct 20 -ct 25";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        // our base file
        File baseOutputFile = this.createTempFile("depthofcoveragenofiltering",".tmp");
        spec.setOutputFileLocation(baseOutputFile);

        // now add the expected files that get generated
        spec.addAuxFile("2f072fd8b41b5ac1108797f89376c797", baseOutputFile);
        spec.addAuxFile("d17ac7cc0b58ba801d2b0727a363d615", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_counts"));
        spec.addAuxFile("c05190c9e6239cdb1cd486edcbc23505", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_proportions"));
        spec.addAuxFile("9cd395f47b329b9dd00ad024fcac9929", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_statistics"));
        spec.addAuxFile("c94a52b4e73a7995319e0b570c80d2f7", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));
        spec.addAuxFile("1970a44efb7ace4e51a37f0bd2dc84d1", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_statistics"));
        spec.addAuxFile("c321c542be25359d2e26d45cbeb6d7ab", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_summary"));
        spec.addAuxFile("9023cc8939777d515cd2895919a99688", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_counts"));
        spec.addAuxFile("3597b69e90742c5dd7c83fbc74d079f3", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_proportions"));
        spec.addAuxFile("7b9d0e93bf5b5313995be7010ef1f528", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("1a6ea3aa759fb154ccc4e171ebca9d02", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_summary"));
        spec.addAuxFile("b492644ff06b4ffb044d5075cd168abf", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_statistics"));
        spec.addAuxFile("77cef87dc4083a7b60b7a7b38b4c0bd8", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_summary"));
        spec.addAuxFile("8e1adbe37b98bb2271ba13932d5c947f", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("761d2f9daf2ebaf43abf65c8fd2fcd05", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("df0ba76e0e6082c0d29fcfd68efc6b77", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("0582b4681dbc02ece2dfe2752dcfd228", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("0685214965bf1863f7ce8de2e38af060", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("7a0cd8a5ebaaa82621fd3b5aed9c32fe", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));

        execute("testBaseOutputNoFiltering",spec);
    }

    @Test
    public void testNoCoverageDueToFiltering() {
        String[] intervals = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/fhs_jhs_30_targts.interval_list"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/FHS_indexed_subset.bam"};

        String cmd = buildRootCmd(hg18Reference,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 5 --maxBaseQuality 4 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        File baseOutputFile = this.createTempFile("depthofcoveragenofiltering",".tmp");
        spec.setOutputFileLocation(baseOutputFile);

        spec.addAuxFile("6ccd7d8970ba98cb95fe41636a070c1c",baseOutputFile);
        spec.addAuxFile("7d87783b3d98b928cac16d383ceca807",createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));

        execute("testNoCoverageDueToFiltering",spec);
    }
}
