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

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ClipReadsWalkersIntegrationTest extends WalkerTest {
    public void testClipper(String name, String args, String md51, String md52) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T ClipReads " +
                        "-I " + validationDataLocation + "clippingReadsTest.bam " +
                        "-o %s " +
                        "-ob %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    final static String Q10ClipOutput = "b29c5bc1cb9006ed9306d826a11d444f";
    @Test public void testQClip0() { testClipper("clipQSum0", "-QT 0", "117a4760b54308f81789c39b1c9de578", "4deca83d80dfa3f093e0dc27d27d1352"); }
    @Test public void testQClip2() { testClipper("clipQSum2", "-QT 2", Q10ClipOutput, "1e123233ebd2f35ac3f41b1b7d2c8199"); }
    @Test public void testQClip10() { testClipper("clipQSum10", "-QT 10", "b29c5bc1cb9006ed9306d826a11d444f", "1e123233ebd2f35ac3f41b1b7d2c8199"); }
    @Test public void testQClip20() { testClipper("clipQSum20", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "b950538d2d8fac1bcee11850c452bd6a"); }
    @Test public void testQClip30() { testClipper("clipQSum30", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "b950538d2d8fac1bcee11850c452bd6a"); }

    @Test public void testClipRange1() { testClipper("clipRange1", "-CT 1-5", "b5acd753226e25b1e088838c1aab9117", "9f70540b795f227668dcf78edcb35c09"); }
    @Test public void testClipRange2() { testClipper("clipRange2", "-CT 1-5,11-15", "be4fcad5b666a5540028b774169cbad7", "a22347a741640fc6df92700e0e8d6f61"); }

    @Test public void testClipSeq() { testClipper("clipSeqX", "-X CCCCC", "db199bd06561c9f2122f6ffb07941fbc", "f49e9e61a44115e2be59330259966f53"); }
    @Test public void testClipSeqFile() { testClipper("clipSeqXF", "-XF " + validationDataLocation + "seqsToClip.fasta", "d011a3152b31822475afbe0281491f8d", "5c977f261442ab6122d5198fa4086e67"); }

    @Test public void testClipMulti() { testClipper("clipSeqMulti", "-QT 10 -CT 1-5 -XF " + validationDataLocation + "seqsToClip.fasta -X CCCCC", "a23187bd9bfb06557f799706d98441de", "38d5f33d198aeee7eebec9feb7b11199"); }

    @Test public void testClipNs() { testClipper("testClipNs", "-QT 10 -CR WRITE_NS", Q10ClipOutput, "1e123233ebd2f35ac3f41b1b7d2c8199"); }
    @Test public void testClipQ0s() { testClipper("testClipQs", "-QT 10 -CR WRITE_Q0S", Q10ClipOutput, "d44cab2e3b70f5492a0f5b59f0b80043"); }
    @Test public void testClipSoft() { testClipper("testClipSoft", "-QT 10 -CR SOFTCLIP_BASES", Q10ClipOutput, "b86374a7e6f59e3dd35781e9e8006702"); }

    @Test
    public void testUseOriginalQuals() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T ClipReads" +
                        " -I " + validationDataLocation + "originalQuals.chr1.1-1K.bam" +
                        " -L chr1:1-1,000" +
                        " -OQ -QT 4 -CR WRITE_Q0S" +
                        " -o %s -ob %s",
                2,
                Arrays.asList("55c01ccc2e84481b22d3632cdb06c8ba", "22db22749f811d30216215e047461621"));
        executeTest("clipOriginalQuals", spec);
    }
}
