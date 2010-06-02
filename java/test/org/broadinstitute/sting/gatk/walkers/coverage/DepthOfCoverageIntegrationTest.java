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

        String cmd = buildRootCmd(hg18,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 0 -dels -baseCounts -pt readgroup -pt sample -pt library --outputFormat csv -ct 10 -ct 15 -ct 20 -ct 25";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        // now add the expected files that get generated
        spec.addAuxFile("494c174ea0773bda98297a0cbdc188eb", baseOutputFile);
        spec.addAuxFile("9df5e7e07efeb34926c94a724714c219", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_counts"));
        spec.addAuxFile("b9a7748e5aec4dc06daed893c901c00d", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_proportions"));
        spec.addAuxFile("9cd395f47b329b9dd00ad024fcac9929", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_statistics"));
        spec.addAuxFile("ca95d2508366d32bf91bf0b0009a023a", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));
        spec.addAuxFile("f6dbd74d32a48abe71ce08d300bce983", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_statistics"));
        spec.addAuxFile("7962a7c09c43ff8b339fa52bce51bfca", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_summary"));
        spec.addAuxFile("b82846df660f0aac8429aec57c2a62d6", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_counts"));
        spec.addAuxFile("d32a8c425fadcc4c048bd8b48d0f61e5", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_proportions"));
        spec.addAuxFile("7b9d0e93bf5b5313995be7010ef1f528", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("755463d88222c81b84f99615e7b4cfd6", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_summary"));
        spec.addAuxFile("e70952f241eebb9b5448f2e7cb288131", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_statistics"));
        spec.addAuxFile("51e4c04dfcb4a20c552ca6f013977fa8", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_summary"));
        spec.addAuxFile("d53431022f7387fe9ac47814ab1fcd88", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("650ee3714da7fbad7832c9d4ad49eb51", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("df0ba76e0e6082c0d29fcfd68efc6b77", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("a52395c883ce8f1a62444d214fe37e88", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("a50011571334f17e950ad3ed1149e350", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("f12a5f97b69718333c4987e3beb98f06", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));
        
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
        spec.addAuxFile("0ee40f3e5091536c14e077b77557083a",createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));

        execute("testNoCoverageDueToFiltering",spec);
    }

    public File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
