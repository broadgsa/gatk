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

package org.broadinstitute.sting.gatk.walkers.coverage;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.util.Arrays;

public class CallableLociWalkerIntegrationTest extends WalkerTest {
    final static String commonArgs = "-R " + b36KGReference + " -T CallableLoci -I " + validationDataLocation + "/NA12878.1kg.p2.chr1_10mb_11_mb.SLX.bam -o %s";

    final static String SUMMARY_MD5 = "cd597a8dae35c226a2cb110b1c9f32d5";

    @Test
    public void testCallableLociWalkerBed() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("c86ac1ef404c11d5e5452e020c8f7ce9", SUMMARY_MD5));
        executeTest("formatBed", spec);
    }

    @Test
    public void testCallableLociWalkerPerBase() {
        String gatk_args = commonArgs + " -format STATE_PER_BASE -L 1:10,000,000-11,000,000 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("d8536a55fe5f6fdb1ee6c9511082fdfd", SUMMARY_MD5));
        executeTest("format_state_per_base", spec);
    }
    
    @Test
    public void testCallableLociWalker2() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-10,000,100 -L 1:10,000,110-10,000,120 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("c671f65712d9575b8b3e1f1dbedc146e", "d287510eac04acf5a56f5cde2cba0e4a"));
        executeTest("formatBed by interval", spec);
    }

    @Test
    public void testCallableLociWalker3() {
        String gatk_args = commonArgs + " -format BED -L 1:10,000,000-11,000,000 -minDepth 10 -maxDepth 100 --minBaseQuality 10 --minMappingQuality 20 -summary %s";
        WalkerTestSpec spec = new WalkerTestSpec(gatk_args, 2,
                Arrays.asList("bc966060184bf4605a31da7fe383464e", "d624eda8f6ed14b9251ebeec73e37867"));
        executeTest("formatBed lots of arguments", spec);
    }
}
