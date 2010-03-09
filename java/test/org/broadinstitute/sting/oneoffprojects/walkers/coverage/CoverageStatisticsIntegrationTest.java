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

        String cmd = buildRootCmd(b36,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) + " -mmq 0 -mbq 0 -dels -baseCounts -both";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());

        // now add the expected files that get generated
        spec.addAuxFile("959937a9b0ace520b4b7d9915d708003", baseOutputFile);
        spec.addAuxFile("aff2349d6dc221c08f6c469379aeaedf", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("6476ed0c54a4307a618aa6d3268b050f", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("50870dad272f03f77befb0075baed1cd", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_locus_statistics"));
        spec.addAuxFile("65318c1e73d98a59cc6f817cde12d3d4", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary_statistics"));
        spec.addAuxFile("ef8c3e2ba3fc0da829e10e2d487c00d2", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("223377e07b35e81a394b75b38d8e72ee", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("096f4ed94020327288ea76245ebd6942", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_summary"));
        spec.addAuxFile("06ed004c86f8b2ad8e64a3b42a0d85c5", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_locus_statistics"));
        spec.addAuxFile("43c160ff9d754744728c142709011993", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_statistics"));
        spec.addAuxFile("a374410efe20609c5c4b87a6da7f4d51", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_summary"));
        
        execute("testBaseOutputNoFiltering",spec);
    }

    @Test
    public void testMedianOverRightHandBin() {
        File base = this.createTempFile("depthofcoveragelowbins",".tmp");
        this.setOutputFileLocation(base);
        String[] intervals = {"1:10,000,000-10,000,800","1:10,250,001-10,250,500","1:10,500,001-10,500,300","1:10,750,001-10,750,400"};
        String[] bams = {"/humgen/gsa-hpprojects/GATK/data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam","/broad/1KG/DCC_merged/freeze5/NA19240.pilot2.454.bam"};

        String cmd = buildRootCmd(b36,new ArrayList<String>(Arrays.asList(bams)),new ArrayList<String>(Arrays.asList(intervals))) +
                " -mmq 0 -mbq 0 -dels -baseCounts -both --start 1 --stop 14 --nBins 13";
        WalkerTestSpec spec = new WalkerTestSpec(cmd,0, new ArrayList<String>());
        spec.addAuxFile("959937a9b0ace520b4b7d9915d708003", base);
        spec.addAuxFile("219d643627eedd696bc476aac96376c2", createTempFileFromBase(base.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("dd0225cf1e0b0bd4289b82fd4939f9fd", createTempFileFromBase(base.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("63575a8a2110507e08d421d44d06b327", createTempFileFromBase(base.getAbsolutePath()+".sample_interval_summary"));

        execute("testMedianOverRHBin",spec);

    }

    public File createTempFileFromBase(String name) {
        File fl = new File(name);
        fl.deleteOnExit();
        return fl;
    }
}
