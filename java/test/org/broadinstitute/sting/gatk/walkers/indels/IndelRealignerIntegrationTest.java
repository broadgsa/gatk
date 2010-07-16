package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {

    @Test
    public void testRealignerLod5() {
        String[] md5s = {"9247b1437ec3b9bc96590524f245220c", "18fca887d1eb7dc300e717ae03b9da62"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod5", spec);
    }

    @Test
    public void testRealignerLod50() {
        String[] md5s = {"9247b1437ec3b9bc96590524f245220c", "9537e4f195ce5840136f60fb61201369"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 50 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod50", spec);
    }

    @Test
    public void testRealignerKnownsOnly() {
        String[] md5s = {"83bc0c9a7d8872b552c6cbd994672c3b", "92bf331b672a03d63c26e4d9d3615f5b"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 1.0 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10076000 -compress 1 -targetIntervals " + validationDataLocation + "NA12878.indels.intervals -B knownIndels,VCF," + validationDataLocation + "NA12878.indels.vcf4 -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe -knownsOnly",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner known indels only", spec);
    }
}