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
                        "-I " + validationDataLocation + "clippingReadsTest.withRG.bam " +
                        "-os %s " +
                        "-o %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    final static String Q10ClipOutput = "b29c5bc1cb9006ed9306d826a11d444f";
    @Test public void testQClip0() { testClipper("clipQSum0", "-QT 0", "117a4760b54308f81789c39b1c9de578", "33e781084379aae538954e30919e8fd3"); }
    @Test public void testQClip2() { testClipper("clipQSum2", "-QT 2", Q10ClipOutput, "57c05b6241db7110148a91fde2d431d0"); }
    @Test public void testQClip10() { testClipper("clipQSum10", "-QT 10", "b29c5bc1cb9006ed9306d826a11d444f", "57c05b6241db7110148a91fde2d431d0"); }
    @Test public void testQClip20() { testClipper("clipQSum20", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "67263a39d5127f2660a5b638ff32056a"); }

    @Test public void testClipRange1() { testClipper("clipRange1", "-CT 1-5", "b5acd753226e25b1e088838c1aab9117", "764846d0592f346a33525af674fd7a10"); }
    @Test public void testClipRange2() { testClipper("clipRange2", "-CT 1-5,11-15", "be4fcad5b666a5540028b774169cbad7", "3061cf742f9e5526a61130128ae761a3"); }

    @Test public void testClipSeq() { testClipper("clipSeqX", "-X CCCCC", "db199bd06561c9f2122f6ffb07941fbc", "b89459f373e40f0b835c1faff2208839"); }
    @Test public void testClipSeqFile() { testClipper("clipSeqXF", "-XF " + validationDataLocation + "seqsToClip.fasta", "d011a3152b31822475afbe0281491f8d", "24e19116ef16a37a6d095ed5c22c2466"); }

    @Test public void testClipMulti() { testClipper("clipSeqMulti", "-QT 10 -CT 1-5 -XF " + validationDataLocation + "seqsToClip.fasta -X CCCCC", "a23187bd9bfb06557f799706d98441de", "ad8d30300cb43d5e300fcc4d2450da8e"); }

    @Test public void testClipNs() { testClipper("testClipNs", "-QT 10 -CR WRITE_NS", Q10ClipOutput, "57c05b6241db7110148a91fde2d431d0"); }
    @Test public void testClipQ0s() { testClipper("testClipQs", "-QT 10 -CR WRITE_Q0S", Q10ClipOutput, "2a1a3153e0942ab355fd8a6e082b30e0"); }
    @Test public void testClipSoft() { testClipper("testClipSoft", "-QT 10 -CR SOFTCLIP_BASES", Q10ClipOutput, "50d43d63d8e39f67a87a6359963c6f52"); }

    @Test
    public void testUseOriginalQuals() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T ClipReads" +
                        " -I " + validationDataLocation + "originalQuals.chr1.1-1K.bam" +
                        " -L chr1:1-1,000" +
                        " -OQ -QT 4 -CR WRITE_Q0S" +
                        " -o %s -os %s",
                2,
                Arrays.asList("22db22749f811d30216215e047461621", "55c01ccc2e84481b22d3632cdb06c8ba"));
        executeTest("clipOriginalQuals", spec);
    }
}
