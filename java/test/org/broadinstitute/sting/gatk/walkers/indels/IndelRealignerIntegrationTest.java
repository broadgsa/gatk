package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {
    @Test
    public void testIntervals() {

        String[] md5lod5 = {"67c3fc25e9d192cc5fbfd48ade0efc84", "86778f92b0fa6aa7c26e651c8c1eb320"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IndelRealigner -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("testLod5", spec1);

        String[] md5lod200 = {"96edef86cea95f312ee8295b38227eb8", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IndelRealigner -LOD 200 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("testLod200", spec2);
    }
}
