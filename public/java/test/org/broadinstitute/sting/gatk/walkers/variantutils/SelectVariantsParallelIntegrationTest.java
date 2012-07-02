package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class SelectVariantsParallelIntegrationTest extends WalkerTest {

    private class ParallelSelectTestProvider extends TestDataProvider {
        final String reference;
        final String args;
        final String md5;
        final int nt;

        private ParallelSelectTestProvider(final String reference, final String args, final String md5, final int nt) {
            super(ParallelSelectTestProvider.class);
            this.reference = reference;
            this.args = args;
            this.md5 = md5;
            this.nt = nt;
        }

        public final String getCmdLine() {
            return "-T SelectVariants -R " + reference + " -o %s -L 1 --no_cmdline_in_header -nt " + nt + " " + args;
        }

        public String toString() {
            return String.format("ParallelSelectVariants nt=%d args=%s", nt, args);
        }
    }

    @DataProvider(name = "ParallelSelectTest")
    public Object[][] makeParallelSelectTestProvider() {
        for ( int nt : Arrays.asList(1, 2, 4) ) {
            { // original MAF test
                String testfile = validationDataLocation + "test.filtered.maf_annotated.vcf";
                String samplesFile = validationDataLocation + "SelectVariants.samples.txt";
                String args = " -sn A -se '[CDH]' -sf " + samplesFile + " -env -ef -select 'DP < 250' --variant " + testfile;
                new ParallelSelectTestProvider(b36KGReference, args, "4386fbb258dcef4437495a37f5a83c53", nt);
            }
            { // new tests on b37 using testdir VCF
                final String testfile = privateTestDir + "NA12878.hg19.example1.vcf";
                final String args = "-select 'DP > 30' -V " + testfile;
                new ParallelSelectTestProvider(b37KGReference, args, "c64b45a14d41b1e5cddbe036b47e7519", nt);
            }
        }

        return ParallelSelectTestProvider.getTests(ParallelSelectTestProvider.class);
    }

    @Test(dataProvider = "ParallelSelectTest")
    public void testParallelSelectTestProvider(final ParallelSelectTestProvider cfg) {
        final WalkerTestSpec spec = new WalkerTestSpec( cfg.getCmdLine(), 1, Arrays.asList(cfg.md5) );
        executeTest(cfg.toString(), spec);
    }
}
