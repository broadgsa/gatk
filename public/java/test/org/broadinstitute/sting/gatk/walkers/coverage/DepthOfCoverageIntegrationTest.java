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
        spec.addAuxFile("19e862f7ed3de97f2569803f766b7433", baseOutputFile);
        spec.addAuxFile("c64cc5636d4880b80b71169ed1832cd7", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_counts"));
        spec.addAuxFile("1a8ba07a60e55f9fdadc89d00b1f3394", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_cumulative_coverage_proportions"));
        spec.addAuxFile("0075cead73a901e3a9d07c5d9c2b75f4", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_statistics"));
        spec.addAuxFile("d757be2f953f893e66eff1ef1f0fff4e", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_interval_summary"));
        spec.addAuxFile("de08996729c774590d6a4954c906fe84", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_statistics"));
        spec.addAuxFile("58ad39b100d1f2af7d119f28ba626bfd", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".library_summary"));
        spec.addAuxFile("0b4ce6059e6587ae5a986afbbcc7d783", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_counts"));
        spec.addAuxFile("adc2b2babcdd72a843878acf2d510ca7", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_cumulative_coverage_proportions"));
        spec.addAuxFile("884281c139241c5db3c9f90e8684d084", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_statistics"));
        spec.addAuxFile("b90636cad74ff4f6b9ff9a596e145bd6", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_interval_summary"));
        spec.addAuxFile("ad540b355ef90c566bebaeabd70026d2", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_statistics"));
        spec.addAuxFile("27fe09a02a5b381e0ed633587c0f4b23", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".read_group_summary"));
        spec.addAuxFile("5fcd53b4bd167b5e6d5f92329cf8678e", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_counts"));
        spec.addAuxFile("7a2a19e54f73a8e07de2f020f1f913dd", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_cumulative_coverage_proportions"));
        spec.addAuxFile("852a079c5e9e93e7daad31fd6a9f4a49", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_statistics"));
        spec.addAuxFile("0828762842103edfaf115ef4e50809c6", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_interval_summary"));
        spec.addAuxFile("5c5aeb28419bba1decb17f8a166777f2", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_statistics"));
        spec.addAuxFile("e5fd6216b3d6a751f3a90677b4e5bf3c", createTempFileFromBase(baseOutputFile.getAbsolutePath()+".sample_summary"));

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
