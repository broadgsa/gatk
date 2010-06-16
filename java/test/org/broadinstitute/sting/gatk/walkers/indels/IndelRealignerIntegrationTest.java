package org.broadinstitute.sting.gatk.walkers.indels;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.io.File;

public class IndelRealignerIntegrationTest extends WalkerTest {
    @Test
    public void testRealigner() {

        String[] md5lod5 = {"96edef86cea95f312ee8295b38227eb8", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec1 = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s",
                 2,
                 Arrays.asList(md5lod5));
        executeTest("test Lod5", spec1);

        String[] md5lod200 = {"96edef86cea95f312ee8295b38227eb8", "d4d8ff567b614729ab8c52bd7d6bef48"};
        WalkerTestSpec spec2 = new WalkerTestSpec(
                "-T IndelRealigner -noPG -LOD 200 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + "NA12878.chrom1.SLX.SRP000032.2009_06.bam -L 1:10023800-10332350 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O %s -snps %s",
                 2,
                 Arrays.asList(md5lod200));
        executeTest("test Lod200", spec2);

        String filename1 = "NA12878.chrom1.SLX.SRP000032.2009_06";
        String filename2 = "low_coverage_CEU.chr1.10k-11k";
        WalkerTestSpec spec3 = new WalkerTestSpec(
                "-T IndelRealigner -nway -noPG -LOD 5 -maxConsensuses 100 -greedy 100 -R " + oneKGLocation + "reference/human_b36_both.fasta -I " + validationDataLocation + filename1 + ".bam -I " + validationDataLocation + filename2 + ".bam -L 1:10023900-10024000 -compress 1 -targetIntervals " + validationDataLocation + "cleaner.test.intervals -O /tmp -snps %s",
                 1,
                 Arrays.asList("bd42a4fa66d7ec7a480c2b94313a78d3"));
        File file1 = new File("/tmp/" + filename1 + ".cleaned.bam");
        file1.deleteOnExit();
        spec3.addAuxFile("1ceae553c8aa20681ed0736d4d2b4541", file1);
        File file2 = new File("/tmp/" + filename2 + ".cleaned.bam");
        file2.deleteOnExit();
        spec3.addAuxFile("ce8ddeae5a5aab836ac1dde9448ccb66", file2);
        executeTest("test NWay", spec3);
    }
}