package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {

    @Test
    public void testRealignerLod5() {
        String[] md5lod5 = {"56f1fb75cae706a5a6278257ea2f2598", "18fca887d1eb7dc300e717ae03b9da62"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("test realigner lod5", spec);
    }

    @Test
    public void testRealignerLod50() {
        String[] md5lod50 = {"56f1fb75cae706a5a6278257ea2f2598", "9537e4f195ce5840136f60fb61201369"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 50 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5lod50));
        executeTest("test realigner lod50", spec);
    }
}