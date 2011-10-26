package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class BAQIntegrationTest extends WalkerTest {
    private final static String baseCommand = "-T PrintReads -R " + b36KGReference +
            " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.allTechs.bam" +
            " -o %s" +
            " -L 1:10,000,000-10,100,000";

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing BAQ with SLX, 454, and SOLID data
    //
    // --------------------------------------------------------------------------------------------------------------
    @Test
    public void testPrintReadsNoBAQ() {
        WalkerTestSpec spec = new WalkerTestSpec( baseCommand +" -baq OFF",  1, Arrays.asList("d97340a2bba2c6320d1ebeb86024a27c"));
        executeTest(String.format("testPrintReadsNoBAQ"), spec);
    }

    @Test
    public void testPrintReadsRecalBAQ() {
        WalkerTestSpec spec = new WalkerTestSpec( baseCommand +" -baq RECALCULATE",  1, Arrays.asList("4ac691bde1ba1301a59857694fda6ae2"));
        executeTest(String.format("testPrintReadsRecalBAQ"), spec);
    }
}
