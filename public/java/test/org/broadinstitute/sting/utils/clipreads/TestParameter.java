package org.broadinstitute.sting.utils.clipreads;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 11/28/11
 * Time: 4:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class TestParameter {
    int inputStart;
    int inputStop;
    int substringStart;
    int substringStop;
    String cigar;

    public TestParameter(int InputStart, int InputStop, int SubstringStart, int SubstringStop, String Cigar) {
        inputStart = InputStart;
        inputStop = InputStop;
        substringStart = SubstringStart;
        substringStop = SubstringStop;
        cigar = Cigar;
    }
}
