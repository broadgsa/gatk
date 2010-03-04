package org.broadinstitute.sting.oneoffprojects.walkers.coverage;

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
public class CoverageStatisticsIntegrationTest extends WalkerTest {

    private boolean RUN_TESTS = true;
    private String root = "-T CoverageStatistics ";
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

        String[] intervals = {"1:10,000,000-10,000,800","1:10,250,001-10,250,500","1:10,500,001-10,500,300","1:10,750,001-10,750,400"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam","/broad/1KG/DCC_merged/freeze5/NA19240.pilot2.454.bam"};

        String cmd = buildRootCmd(b36,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 0 -dels -baseCounts";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        // now add the expected files that get generated
        spec.addAuxFile("cb87d6069ac60c73f047efc6d9386619", baseOutputFile);
        spec.addAuxFile("aff2349d6dc221c08f6c469379aeaedf", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".interval_statistics"));
        spec.addAuxFile("6476ed0c54a4307a618aa6d3268b050f", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".interval_summary"));
        spec.addAuxFile("c744a298b7541f3f823e6937e9a0bc67", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".locus_statistics"));
        spec.addAuxFile("65318c1e73d98a59cc6f817cde12d3d4", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".summary_statistics"));
        spec.addAuxFile("9fc19f773a7ddfbb473d124e675a3d94", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        execute("testBaseOutputNoFiltering",spec);
    }

    public File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
