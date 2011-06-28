package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class WalkerTestIntegrationTest extends WalkerTest {

    public void testBadMD5(String md5) {
        WalkerTestSpec spec = new WalkerTestSpec("FAIL", Arrays.asList(md5));
        executeTest("", spec);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testNullMD5() {
        testBadMD5(null);
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testBadLengthMD5() {
        testBadMD5("asdfasdfa");
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testSpacesMD5() {
        testBadMD5("1de8e943fbf55246ebd19efa32f22a58 ");
    }

    @Test(expectedExceptions = RuntimeException.class)
    public void testBadCharMD5() {
        testBadMD5("1de8e943fbf55246ebd19efa32f22a5_");
    }
}
