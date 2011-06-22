package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.HashMap;

public class PrintReadsIntegrationTest extends WalkerTest {
    private static class PRTest {
        final static String REF = hg18Reference;
        final static String BAM = validationDataLocation + "HiSeq.1mb.bam";
        String args;
        String md5;

        private PRTest(String args, String md5) {
            this.args = args;
            this.md5 = md5;
        }
    }

    @DataProvider(name = "PRTest")
    public Object[][] createData1() {
        return new Object[][]{
                {new PRTest("",  "dc8e5451dd29757c336013146010f73a")},
                {new PRTest(" -compress 0",  "fde82269c78c9e91e57286433531b4af")},
                {new PRTest(" -simplifyBAM",  "0531717b32a7e21c0de70b1526b0751f")},
                {new PRTest(" -n 10",  "cdc4ddf9ee1d2ecf37168da8ef23c270")} };
    }

    @Test(dataProvider = "PRTest")
    public void testPrintReads(PRTest params) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads -R " + params.REF +
                        " -I " + params.BAM +
                        params.args +
                        " -o %s",
                Arrays.asList(params.md5));
        executeTest("testPrintReads-"+params.args, spec).getFirst();
    }

    @Test
    public void testPrintReadsReadAllInsertion() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T PrintReads -R " + b37KGReference +
                        " -I " + validationDataLocation + "oneReadAllInsertion.bam" +
                        " -o %s",
                Arrays.asList("6caec4f8a25befb6aba562955401af93"));
        executeTest("testPrintReads-oneReadAllInsertion", spec);
    }
}

