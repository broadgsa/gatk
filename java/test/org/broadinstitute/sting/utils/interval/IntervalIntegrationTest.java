/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.interval;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;

/**
 * Test the GATK core interval parsing mechanism.
 */
public class IntervalIntegrationTest extends WalkerTest {
    @Test
    public void testAllImplicitIntervalParsing() {
        String md5 = "7821db9e14d4f8e07029ff1959cd5a99";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T CountLoci" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -o %s",
                        1, // just one output file
                        Arrays.asList(md5));
        executeTest("testAllIntervalsImplicit",spec);
    }

    @Test
    public void testAllExplicitIntervalParsing() {
        String md5 = "7821db9e14d4f8e07029ff1959cd5a99";
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T CountLoci" +
                        " -I " + validationDataLocation + "OV-0930.normal.chunk.bam" +
                        " -R " + hg18Reference +
                        " -L all" +
                        " -o %s",
                        1, // just one output file
                        Arrays.asList(md5));                        
        executeTest("testAllIntervalsExplicit",spec);
    }
}
