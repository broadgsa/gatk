package org.broadinstitute.sting.alignment;

import org.testng.annotations.Test;
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
        String md5 = "a2bdf907b18114a86ca47f9fc23791bf";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + GATKDataLocation + "bwa/human_b36_both.fasta" +
                        " -T Align" +
                        " -I " + validationDataLocation + "NA12878_Pilot1_20.trimmed.unmapped.bam" +
                        " -o %s",
                1, // just one output file
                Arrays.asList(md5));
        executeTest("testBasicAlignment", spec);
    }
}
