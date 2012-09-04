package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 8/31/12
 * Time: 11:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class InvalidArgumentIntegrationTest extends WalkerTest {
    private static final String callsB36  = BaseTest.validationDataLocation + "lowpass.N3.chr1.raw.vcf";

    private WalkerTest.WalkerTestSpec baseTest(String flag, String arg, Class exeption) {
        return new WalkerTest.WalkerTestSpec("-T VariantsToTable -M 10 --variant:vcf "
                + callsB36 + " -F POS,CHROM -R "
                + b36KGReference +  " -o %s " + flag + " " + arg,
                1, exeption);

    }

    @Test
    public void testUnknownReadFilter() {
        executeTest("UnknownReadFilter",baseTest("-rf","TestUnknownReadFilter", UserException.MalformedReadFilterException.class));
    }

    @Test
    public void testMalformedWalkerArgs() {
        executeTest("MalformedWalkerArgs",
                new WalkerTest.WalkerTestSpec("-T UnknownWalkerName -M 10 --variant:vcf "
                + callsB36 + " -F POS,CHROM -R "
                + b36KGReference +  " -o %s ",
                1, UserException.MalformedWalkerArgumentsException.class));
    }
}
