package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.List;
import java.io.File;

public class DepthOfCoverageIntegrationTest extends WalkerTest {
    private static String root = "-L 1:10,164,500-10,164,520 -R /broad/1KG/reference/human_b36_both.fasta -T DepthOfCoverage -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam";
    static HashMap<String, String> expectations = new HashMap<String, String>();
    static {
        expectations.put("-minMAPQ 1", "59c6071105a598e19f460640c35768c6");
        expectations.put("-minMAPQ 100", "e997fb5d61eaec21518722b0de90af20");
        expectations.put("-minDepth 8", "3e50afef0e751119cd27c324bdfae544");
        expectations.put("-minDepth 10", "d4c336d9e748347e1082bbc92d2489a3");
        expectations.put("-bySample", "160ffa185dbfa8b0d2dc57f60f5b1e48");
        expectations.put("-byRG", "dd3b4d040df7325dad4760ac6fa5252d");
        expectations.put("-minMAPQ 1 -bySample -byRG", "bd2a07ef548b86e82ac6cce534225612");
    }

    @Test
    public void testDepthOfCoverage1() {

        for ( Map.Entry<String, String> entry : expectations.entrySet() ) {
            String extraArgs = entry.getKey();
            String md5 = entry.getValue();

            WalkerTestSpec spec = new WalkerTestSpec( root + " " + extraArgs + " -o %s",
                    1, // just one output file
                    Arrays.asList(md5));
            executeTest("testDepthOfCoverage1", spec);
        }
    }

    @Test
    public void testDepthOfCoverage454() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T DepthOfCoverage -R /broad/1KG/reference/human_b36_both.fasta -I /humgen/gsa-scr1/GATK_Data/Validation_Data/NA12873.454.SRP000031.2009_06.chr1.10_20mb.bam -L 1:10,001,890-10,001,895 -o %s",
                1, // just one output file
                Arrays.asList("51203ba5ab928449cd01363af0b91510"));
        executeTest("testDepthOfCoverage454", spec);
    }    
}