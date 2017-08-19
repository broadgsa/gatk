/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.readutils;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ClipReadsWalkersIntegrationTest extends WalkerTest {
    public void testClipper(String name, String args, String md51, String md52) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T ClipReads " +
                        "-I " + privateTestDir + "clippingReadsTest.withRG.bam " +
                        "-os %s " +
                        "-o %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    final static String Q10ClipOutput = "b29c5bc1cb9006ed9306d826a11d444f";
    @Test public void testQClip0() { testClipper("clipQSum0", "-QT 0", "117a4760b54308f81789c39b1c9de578", "2e617e0f0b4a1bf45dbb666da0a6a16c"); }
    @Test public void testQClip2() { testClipper("clipQSum2", "-QT 2", Q10ClipOutput, "1b0d2317dfd1fd5a4d14eea3f560ac17"); }
    @Test public void testQClip10() { testClipper("clipQSum10", "-QT 10", "b29c5bc1cb9006ed9306d826a11d444f", "1b0d2317dfd1fd5a4d14eea3f560ac17"); }
    @Test public void testQClip20() { testClipper("clipQSum20", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "3a9a10966cce58baf9e9400a2c5f99ca"); }

    @Test public void testClipRange1() { testClipper("clipRange1", "-CT 1-5", "b5acd753226e25b1e088838c1aab9117", "701f03388ffaeb0f92be1e719e2477df"); }
    @Test public void testClipRange2() { testClipper("clipRange2", "-CT 1-5,11-15", "be4fcad5b666a5540028b774169cbad7", "7dc5da735c437b6523e0823a1ed0e981"); }

    @Test public void testClipSeq() { testClipper("clipSeqX", "-X CCCCC", "db199bd06561c9f2122f6ffb07941fbc", "daa41acd8c9f1a75c8107d0fd0db8cf0"); }
    @Test public void testClipSeqFile() { testClipper("clipSeqXF", "-XF " + privateTestDir + "seqsToClip.fasta", "22cd33dbd06cc425ce3626360367bfca", "b1812de69c967059362139a269044362"); }

    @Test public void testClipMulti() { testClipper("clipSeqMulti", "-QT 10 -CT 1-5 -XF " + privateTestDir + "seqsToClip.fasta -X CCCCC", "e4a88276a6f76ccc2b05a06c76305a1c", "9f999ed70684fdfb8dba019ff71741fe"); }

    @Test public void testClipNs() { testClipper("testClipNs", "-QT 10 -CR WRITE_NS", Q10ClipOutput, "1b0d2317dfd1fd5a4d14eea3f560ac17"); }
    @Test public void testClipQ0s() { testClipper("testClipQs", "-QT 10 -CR WRITE_Q0S", Q10ClipOutput, "b7ce95f277d5f66b7f5a60ffbf9692cf"); }
    @Test public void testClipSoft() { testClipper("testClipSoft", "-QT 10 -CR SOFTCLIP_BASES", Q10ClipOutput, "55c9d86c18c5a59a98429e8a5076e871"); }

    @Test
    public void testUseOriginalQuals() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + hg18Reference +
                        " -T ClipReads" +
                        " -I " + privateTestDir + "originalQuals.chr1.1-1K.bam" +
                        " -L chr1:1-1,000" +
                        " -OQ -QT 4 -CR WRITE_Q0S" +
                        " -o %s -os %s",
                2,
                Arrays.asList("0e4c116fa07c3d8190eefe6968a0b819", "55c01ccc2e84481b22d3632cdb06c8ba"));
        executeTest("clipOriginalQuals", spec);
    }
}
