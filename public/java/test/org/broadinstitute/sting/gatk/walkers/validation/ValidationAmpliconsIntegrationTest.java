package org.broadinstitute.sting.gatk.walkers.validation;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: 7/19/11
 * Time: 7:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class ValidationAmpliconsIntegrationTest extends WalkerTest {

    @Test(enabled=true)
    public void testWikiExample() {
        String siteVCF = validationDataLocation + "sites_to_validate.vcf";
        String maskVCF = validationDataLocation + "amplicon_mask_sites.vcf";
        String intervalTable = validationDataLocation + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("27f9450afa132888a8994167f0035fd7"));
        executeTest("Test probes", spec);
    }

    @Test(enabled=true)
    public void testWikiExampleNoBWA() {
        String siteVCF = validationDataLocation + "sites_to_validate.vcf";
        String maskVCF = validationDataLocation + "amplicon_mask_sites.vcf";
        String intervalTable = validationDataLocation + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30 --doNotUseBWA";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("f2611ff1d9cd5bedaad003251fed8bc1"));
        executeTest("Test probes", spec);
    }

    @Test(enabled=true)
    public void testWikiExampleMonoFilter() {
        String siteVCF = validationDataLocation + "sites_to_validate.vcf";
        String maskVCF = validationDataLocation + "amplicon_mask_sites.vcf";
        String intervalTable = validationDataLocation + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30 --filterMonomorphic";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("77b3f30e38fedad812125bdf6cf3255f"));
        executeTest("Test probes", spec);
    }

}
