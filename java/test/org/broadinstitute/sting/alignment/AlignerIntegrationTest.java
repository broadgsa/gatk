package org.broadinstitute.sting.alignment;

import org.junit.Test;
import org.broadinstitute.sting.WalkerTest;

import java.util.Arrays;

/**
 * Integration tests for the aligner.
 *
 * @author mhanna
 * @version 0.1
 */
public class AlignerIntegrationTest extends WalkerTest {
    @Test
    public void testBasicAlignment() {
        String md5 = "34eb4323742999d6d250a0aaa803c6d5";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R /humgen/gsa-scr1/GATK_Data/bwa/human_b36_both.fasta" +
                        " -T Align" +
                        " -I " + validationDataLocation + "NA12878_Pilot1_20.trimmed.unmapped.bam" +
                        " -ob %s",
                1, // just one output file
                Arrays.asList(md5));
        //executeTest("testBasicAlignment", spec);
        System.err.println("*********************************************************");
        System.err.println("FIX THIS TEST SE TEAM");
        System.err.println("*********************************************************");
                
    }
}
