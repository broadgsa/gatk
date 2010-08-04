/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
 
public class BeagleIntegrationTest extends WalkerTest {

    private static final String beagleValidationDataLocation = validationDataLocation + "/Beagle/";
    @Test
    public void testBeagleOutput() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T BeagleOutputToVCF -R " + hg19Reference + " " +
                        "-B variant,VCF," + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                        "-B beagleR2,BEAGLE," + beagleValidationDataLocation + "inttestbgl.r2 " +
                        "-B beagleProbs,BEAGLE," + beagleValidationDataLocation + "inttestbgl.gprobs " +
                        "-B beaglePhased,BEAGLE," + beagleValidationDataLocation + "inttestbgl.phased " +
                        "-output %s", 1, Arrays.asList("e7b9aac20246f26ffcc599850f6cb2a0"));
        executeTest("test BeagleOutputToVCF", spec);
    }
   
    @Test
    public void testBeagleInput() {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T ProduceBeagleInput -R " + hg19Reference + " " +
                        "-B variant,VCF," + beagleValidationDataLocation + "inttestbgl.input.vcf " +
                         "-beagle %s", 1, Arrays.asList("c0d30e5dbe903874f8422a0e63a5118e"));
        executeTest("test BeagleInput", spec);
    }

}