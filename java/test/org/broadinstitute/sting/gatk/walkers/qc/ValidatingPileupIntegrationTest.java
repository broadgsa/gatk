package org.broadinstitute.sting.gatk.walkers.qc;

import org.junit.Test;
import org.broadinstitute.sting.WalkerTest;

import java.util.Collections;

/**
 * Run validating pileup across a set of core data as proof of the integrity of the GATK core.
 *
 * @author mhanna
 * @version 0.1
 */
public class ValidatingPileupIntegrationTest extends WalkerTest {
    @Test
    public void testEcoliThreaded() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T ValidatingPileup" +
                " -I " + validationDataLocation + "MV1994.selected.bam" +
                " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                " -B:pileup,SAMPileup "+ validationDataLocation + "MV1994.selected.pileup" +
                " -S SILENT -nt 8 -U ALLOW_SEQ_DICT_INCOMPATIBILITY",0, Collections.<String>emptyList());
        executeTest("testEcoliThreaded",spec);
    }
}
