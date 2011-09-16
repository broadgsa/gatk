/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.WalkerTest;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.Test;

/**
 *
 */
public class EngineFeaturesIntegrationTest extends WalkerTest {
    private void testBadRODBindingInput(String type, String name, Class c) {
        WalkerTestSpec spec = new WalkerTestSpec("-T SelectVariants -L 1:1 --variant:variant," + type + " "
                + b37dbSNP132 + " -R " + b37KGReference + " -o %s",
                1, c);
        executeTest(name, spec);
    }

    @Test() private void testBadRODBindingInputType1() {
        testBadRODBindingInput("beagle", "BEAGLE input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputType2() {
        testBadRODBindingInput("vcf3", "VCF3 input to VCF expecting walker", UserException.class);
    }

    @Test() private void testBadRODBindingInputType3() {
        testBadRODBindingInput("bed", "Bed input to VCF expecting walker", UserException.BadArgumentValue.class);
    }

    @Test() private void testBadRODBindingInputTypeUnknownType() {
        testBadRODBindingInput("bedXXX", "Unknown input to VCF expecting walker", UserException.UnknownTribbleType.class);
    }
}