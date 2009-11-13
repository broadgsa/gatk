package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IntervalCleanerIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5lod5 = {"163d2f1b04741986151636d3c7e04c94", "4aa3637e86822c95af3e2c9b414530c3"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("testLod5", spec1);

        String[] md5lod200 = {"1d89ee2af03df79eb5de494c77767221", "e39aa718c9810364ebe30964d878d5ff"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 200 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("testLod200", spec2);

        String[] md5cleanedOnly = {"710f9114ec496c73f1c3782b5cb09757", "4aa3637e86822c95af3e2c9b414530c3"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -cleanedOnly -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5cleanedOnly));
        executeTest("testCleanedOnly", spec3);
    }
}