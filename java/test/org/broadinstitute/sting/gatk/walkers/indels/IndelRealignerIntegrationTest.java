package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

public class IndelRealignerIntegrationTest extends WalkerTest {

    @Test
    public void testRealignerLod5() {
        String[] md5s = {"a377de4e2eb4df8ef79590e4131afe35", "c4ef635f2597b12b93a73199f07e509b"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -o %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod5", spec);
    }

    @Test
    public void testRealignerLod50() {
        String[] md5s = {"a377de4e2eb4df8ef79590e4131afe35", "3735a510513b6fa4161d92155e026283"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 50 -maxConsensuses 100 -greedy 100 -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10030000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -o %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner lod50", spec);
    }

    @Test
    public void testRealignerKnownsOnly() {
        String[] md5s = {"654cb0c845c3f25af6a3f8911ac06a73", "74652bd8240291293ec921f8ecfa1622"};
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 1.0 -R " + b36KGReference + " -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023000-10076000 -compress 1 -targetIntervals " + validationDataLocation + "NA12878.indels.intervals -B:knownIndels,VCF " + validationDataLocation + "NA12878.indels.vcf4 -o %s -stats %s --sortInCoordinateOrderEvenThoughItIsHighlyUnsafe -knownsOnly",
                 2,
                 Arrays.asList(md5s));
        executeTest("test realigner known indels only", spec);
    }
}