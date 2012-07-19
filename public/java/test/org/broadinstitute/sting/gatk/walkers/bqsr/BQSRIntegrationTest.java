package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author ebanks
 * @since 7/16/12
 */
public class BQSRIntegrationTest extends WalkerTest {

    private static class BQSRTest {
        final String reference;
        final String interval;
        final String bam;
        final String args;
        final String md5;

        private BQSRTest(String reference, String bam, String interval, String args, String md5) {
            this.reference = reference;
            this.bam = bam;
            this.interval = interval;
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    @DataProvider(name = "BQSRTest")
    public Object[][] createBQSRTestData() {
        String HiSeqBam = publicTestDir + "HiSeq.1mb.1RG.bam";
        String HiSeqInterval = "chr1:10,000,000-10,100,000";
        return new Object[][]{
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, "", "087dbc3e3f0cee6b891aecad2d0faf25")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov ContextCovariate", "8a69591f728b3a2cdd79ff26bbebcc26")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --no_standard_covs -cov CycleCovariate", "73d649bce0b69f56452de8c7e0a8686d")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --indels_context_size 4", "d9512cebf54ea120539059976b33d1ca")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --low_quality_tail 5", "f61a8df03aae8c4273acf2b72497f154")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --quantizing_levels 6", "7c2ce84e521d8f19fe5061b4e40317f7")},
                {new BQSRTest(hg18Reference, HiSeqBam, HiSeqInterval, " --mismatches_context_size 4", "66a0caad65ab41d9013e812617e67370")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", "", "6f5e9836147b488a7a204cffa5ecd21e")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-10,200,000", "", "444fdfca7835e6a3714445f7e60abcaf")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12873.454.SRP000031.2009_06.chr1.10_20mb.1RG.bam", "1:10,000,000-10,200,000", "", "e0bfaf38f45142d45c8fe0ae05d0d4e0")},
                {new BQSRTest(b36KGReference, validationDataLocation + "originalQuals.1kg.chr1.1-1K.1RG.bam", "1:1-1,000", " -OQ", "5b30760bab51b4d1fc02097d4eacefa4")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA19240.chr1.BFAST.SOLID.bam", "1:10,000,000-20,000,000", " --solid_recal_mode REMOVE_REF_BIAS", "742fd8edfa36ab9023ceeaac40c4e215")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:anyNameABCD,VCF " + privateTestDir + "vcfexample3.vcf", "6f5e9836147b488a7a204cffa5ecd21e")},
                {new BQSRTest(b36KGReference, validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.1Mb.1RG.bam", "1:10,000,000-10,200,000", " -knownSites:bed " + validationDataLocation + "bqsrKnownTest.bed", "22dd42897cf20852712c6e8f63195443")},
        };
    }

    @Test(dataProvider = "BQSRTest")
    public void testBQSR(BQSRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T BaseQualityScoreRecalibrator" +
                        " -R " + params.reference +
                        " -I " + params.bam +
                        " -L " + params.interval +
                        params.args +
                        " --no_plots" +
                        " -knownSites " + (params.reference.equals(b36KGReference) ? b36dbSNP129 : hg18dbSNP132) +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testBQSR-"+params.args, spec).getFirst();
    }

    @Test
    public void testBQSRFailWithoutDBSNP() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + b36KGReference +
                        " -T BaseQualityScoreRecalibrator" +
                        " -I " + validationDataLocation + "NA12892.SLX.SRP000031.2009_06.selected.bam" +
                        " -L 1:10,000,000-10,200,000" +
                        " --no_plots" +
                        " -o %s",
                1, // just one output file
                UserException.CommandLineException.class);
        executeTest("testBQSRFailWithoutDBSNP", spec);
    }

    private static class PRTest {
        final String args;
        final String md5;

        private PRTest(String args, String md5) {
            this.args = args;
            this.md5 = md5;
        }

        @Override
        public String toString() {
            return String.format("PrintReads(args='%s')", args);
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createPRTestData() {
        return new Object[][]{
                {new PRTest("", "d2d6ed8667cdba7e56f5db97d6262676")},
                {new PRTest(" -qq -1", "b7053d3d67aba6d8892f0a60f0ded338")},
                {new PRTest(" -qq 6", "bfbf0855185b2b70aa35237fb71e4487")},
                {new PRTest(" -DIQ", "66aa65223f192ee39c1773aa187fd493")}
        };
    }

    @Test(dataProvider = "PRTest")
    public void testPR(PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads" +
                        " -R " + hg18Reference +
                        " -I " + publicTestDir + "HiSeq.1mb.1RG.bam" +
                        " -BQSR " + publicTestDir + "HiSeq.1mb.1RG.table" +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }
}
