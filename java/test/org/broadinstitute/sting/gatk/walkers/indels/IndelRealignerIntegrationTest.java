package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {
    @Test
    public void testRealigner() {

        String[] md5lod5 = {"d9cbff4832fc3ee7a7ad1c58cc891bdd", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("test Lod5", spec1);

        String[] md5lod200 = {"d9cbff4832fc3ee7a7ad1c58cc891bdd", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 200 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("test Lod200", spec2);
    }
}