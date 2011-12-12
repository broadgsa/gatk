package org.broadinstitute.sting.utils.clipreads;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 11/29/11
 * Time: 4:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class CigarStringTestPair {
    public String toTest;
    public String expected;

    public CigarStringTestPair(String ToTest, String Expected) {
        this.toTest = ToTest;
        this.expected = Expected;
    }
}
