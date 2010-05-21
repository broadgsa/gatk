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
 * @Date Feb 25, 2010
 */
public class DepthOfCoverageIntegrationTest extends WalkerTest {

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

    private void execute(String name, WalkerTestSpec spec) {
        if ( RUN_TESTS ) {
            executeTest(name,spec);
        }
    }

    @Test
    public void testBaseOutputNoFiltering() {
        // our base file
        File baseOutputFile = this.createTempFile("depthofcoveragenofiltering",".tmp");
        this.setOutputFileLocation(baseOutputFile);

        String[] intervals = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/fhs_jhs_30_targts.interval_list"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/FHS_indexed_subset.bam"};

        String cmd = buildRootCmd(hg18,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 0 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        // now add the expected files that get generated
        spec.addAuxFile("fc742e346be2344557cf8c039f467508", baseOutputFile);
        spec.addAuxFile("e58b701b01ec0dbe75c146295434ba3b", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_counts"));
        spec.addAuxFile("b9a7748e5aec4dc06daed893c901c00d", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_proportions"));
        spec.addAuxFile("848e556ec7e03e9b0398d189d7cbb4ad", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_statistics"));
        spec.addAuxFile("e6fc8f7a5fcc440e21de5891f3403d5d", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));
        spec.addAuxFile("cac8e7a688d9bbe781232c61091d3237", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_statistics"));
        spec.addAuxFile("e452dfb5581a7aafaf2122c5ae497f1b", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_summary"));
        spec.addAuxFile("38fb89e1bb52d0342f97f72e86956b79", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_counts"));
        spec.addAuxFile("f9f2941ee39577ac2f80668e7f6b3d4b", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_proportions"));
        spec.addAuxFile("7b9d0e93bf5b5313995be7010ef1f528", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("fd29fe0c14351e934a6fef9df1f38f7b", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_summary"));
        spec.addAuxFile("cc7ee5075a932dba486e78824ca34202", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_statistics"));
        spec.addAuxFile("e1653480daa2d96f7c584ed4cd20c147", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_summary"));
        spec.addAuxFile("529353375d23c529228b38119c51e269", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("650ee3714da7fbad7832c9d4ad49eb51", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("925cc5b49286e0222bce6251d1baafc7", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("d9e812398d778f28ed12d7f3d18628e2", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("f3315551081331bc322c53b11412d707", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("3a059ad82d945643dc4e03239c4041f5", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));
        
        execute("testBaseOutputNoFiltering",spec);
    }

    @Test
    public void testNoCoverageDueToFiltering() {
        File baseOutputFile = this.createTempFile("depthofcoveragenofiltering",".tmp");
        this.setOutputFileLocation(baseOutputFile);

        String[] intervals = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/fhs_jhs_30_targts.interval_list"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/FHS_indexed_subset.bam"};

        String cmd = buildRootCmd(hg18,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 5 --maxBaseQuality 4 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        spec.addAuxFile("d570c27d82a80ebd2852e9d34aff4e87",baseOutputFile);
        spec.addAuxFile("00107eea991f7379771b29dea0c859cb",createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));

        execute("testNoCoverageDueToFiltering",spec);
    }

    public File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
