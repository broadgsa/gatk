package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IntervalCleanerIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5lod5 = {"042c1c2649a51a260bc204ec5f256c1a", "460631e8d98644dcf53b3045ca40f02a"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("testLod5", spec1);

        String[] md5lod200 = {"1d89ee2af03df79eb5de494c77767221", "6137bf0c25c7972b07b0d3fc6979cf5b"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 200 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("testLod200", spec2);

        String[] md5cleanedOnly = {"168c8d02fceb107477381b93d189fe1f", "460631e8d98644dcf53b3045ca40f02a"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -cleanedOnly -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5cleanedOnly));
        executeTest("testCleanedOnly", spec3);
    }
}