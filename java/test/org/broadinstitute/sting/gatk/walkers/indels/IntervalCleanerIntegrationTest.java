package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IntervalCleanerIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5lod5 = {"bfe1c76bf352b22f79c9b7242197a126", "4aa3637e86822c95af3e2c9b414530c3"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IntervalCleaner -stats blah1.stats -LOD 5 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("testLod5", spec1);

        String[] md5lod200 = {"4481e4d24d61a3e438323b368ad0eee7", "e39aa718c9810364ebe30964d878d5ff"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 200 -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("testLod200", spec2);

        String[] md5cleanedOnly = {"c8da42cc2298e05ab4f2bc9e5d7445e1", "4aa3637e86822c95af3e2c9b414530c3"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -cleanedOnly -maxConsensuses 100 -greedy 100 -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.chrom1.SLX.SRP000032.2009_06.bam -L /humgen/gsa-scr1/GATK_Data/Validation_Data/cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5cleanedOnly));
        executeTest("testCleanedOnly", spec3);
    }
}