package org.broadinstitute.sting.playground.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;

/**
 * The pile-up tests, that test any changes to the underlying ROD system
 */
public class RodSystemValidationIntegrationTest extends WalkerTest {

    public static String baseTestString() {
        return "-T RodSystemValidation -o %s -R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta";
    }


    @Test
    public void testSimpleGeliPileup() {
        WalkerTestSpec spec = new WalkerTestSpec(
                baseTestString() + " -B eval,Variants," + validationDataLocation + "ROD_validation/chr1.geli", 1,
                Arrays.asList("536567b13ea4b8786badd96c879df245"));
        executeTest("testVCFSelect1", spec);
    }
}
