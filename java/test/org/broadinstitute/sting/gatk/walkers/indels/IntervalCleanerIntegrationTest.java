package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IntervalCleanerIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5lod5 = {"3c4339cb2a258d1c67bf4d930fe0055a", "86778f92b0fa6aa7c26e651c8c1eb320"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L " + validationDataLocation + "cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("testLod5", spec1);

        String[] md5lod200 = {"32401cef2134d973ff0037df27f1dcca", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 200 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L " + validationDataLocation + "cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("testLod200", spec2);

        String[] md5cleanedOnly = {"6b03a934c65a4964d29c147de2f89144", "86778f92b0fa6aa7c26e651c8c1eb320"};
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T IntervalCleaner -LOD 5 -cleanedOnly -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L " + validationDataLocation + "cleaner.test.intervals --OutputCleaned %s -snps %s",
                 2,
                 Arrays.asList(md5cleanedOnly));
        executeTest("testCleanedOnly", spec3);
    }
}
