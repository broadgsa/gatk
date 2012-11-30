package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.classloader.GATKLiteUtils;
import org.testng.SkipException;
import org.testng.annotations.Test;

import java.util.Arrays;

// ********************************************************************************** //
// Note that this class also serves as an integration test for the VariantAnnotator!  //
// ********************************************************************************** //

public class UnifiedGenotyperLiteIntegrationTest extends WalkerTest {

    private final static String baseCommand = "-T UnifiedGenotyper -R " + b36KGReference + " --no_cmdline_in_header -glm BOTH -minIndelFrac 0.0 --dbsnp " + b36dbSNP129;

    // --------------------------------------------------------------------------------------------------------------
    //
    // testing contamination down-sampling gets ignored
    //
    // --------------------------------------------------------------------------------------------------------------

    @Test
    public void testContaminationDownsampling() {
        if ( !GATKLiteUtils.isGATKLite() )
            throw new SkipException("Only want to test for GATK lite");

        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                baseCommand + " -I " + validationDataLocation + "NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s -L 1:10,000,000-10,010,000", 1,
                Arrays.asList("9addd225a985178339a0c49dc5fdc220"));
        executeTest("test contamination_percentage_to_filter gets ignored", spec);
    }

}
