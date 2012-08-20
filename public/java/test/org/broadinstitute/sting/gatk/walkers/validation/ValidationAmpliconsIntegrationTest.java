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
        String maskVCF = privateTestDir + "amplicon_mask_sites.vcf";
        String intervalTable = privateTestDir + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("240d99b58f73985fb114abe9044c0271"));
        executeTest("Test probes", spec);
    }

    @Test(enabled=true)
    public void testWikiExampleNoBWA() {
        String siteVCF = privateTestDir + "sites_to_validate.vcf";
        String maskVCF = privateTestDir + "amplicon_mask_sites.vcf";
        String intervalTable = privateTestDir + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30 --doNotUseBWA";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("6e7789445e29d91979a21e78d3d53295"));
        executeTest("Test probes", spec);
    }

    @Test(enabled=true)
    public void testWikiExampleMonoFilter() {
        String siteVCF = privateTestDir + "sites_to_validate.vcf";
        String maskVCF = privateTestDir + "amplicon_mask_sites.vcf";
        String intervalTable = privateTestDir + "amplicon_interval_table1.table";
        String testArgs = "-R " + b37KGReference + " -T ValidationAmplicons --ValidateAlleles:VCF "+siteVCF+" -o %s";
        testArgs += " --ProbeIntervals:table "+intervalTable+" -L:table "+intervalTable+" --MaskAlleles:VCF "+maskVCF;
        testArgs += " --virtualPrimerSize 30 --filterMonomorphic";
        WalkerTestSpec spec = new WalkerTestSpec(testArgs, 1,
                Arrays.asList("18d7236208db603e143b40db06ef2aca"));
        executeTest("Test probes", spec);
    }

}
