package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {

    @Test
    public void testRealignerLod5() {
        String[] md5s = {"ce789182974609da27f3eadbbbcd98c4", "18fca887d1eb7dc300e717ae03b9da62"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod5", spec);
    }

    @Test
    public void testRealignerLod50() {
        String[] md5s = {"ce789182974609da27f3eadbbbcd98c4", "9537e4f195ce5840136f60fb61201369"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 50 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod50", spec);
    }

    @Test
    public void testRealignerKnownsOnly() {
        String[] md5s = {"f50b9bddd78ed6fb820d2cbdfbff8300", "1091436c40d5ba557d85662999cc0c1d"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 1.0 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10076000 -compress 1 -targetIntervals " + validationDataLocation + "NA12878.indels.intervals -B knownIndels,VCF," + validationDataLocation + "NA12878.indels.vcf -O %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe -knownsOnly",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner known indels only", spec);
    }
}