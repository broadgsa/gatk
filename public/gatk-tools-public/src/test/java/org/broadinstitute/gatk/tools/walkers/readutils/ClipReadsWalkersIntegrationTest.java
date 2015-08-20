/*
* Copyright 2012-2015 Broad Institute, Inc.
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
    @Test public void testQClip0() { testClipper("clipQSum0", "-QT 0", "117a4760b54308f81789c39b1c9de578", "bcf0d1e13537f764f006ef6d9b401ea7"); }
    @Test public void testQClip2() { testClipper("clipQSum2", "-QT 2", Q10ClipOutput, "27847d330b962e60650df23b6efc8c3c"); }
    @Test public void testQClip10() { testClipper("clipQSum10", "-QT 10", "b29c5bc1cb9006ed9306d826a11d444f", "27847d330b962e60650df23b6efc8c3c"); }
    @Test public void testQClip20() { testClipper("clipQSum20", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "f89ec5439e88f5a75433150da0069034"); }

    @Test public void testClipRange1() { testClipper("clipRange1", "-CT 1-5", "b5acd753226e25b1e088838c1aab9117", "987007f6e430cad4cb4a8d1cc1f45d91"); }
    @Test public void testClipRange2() { testClipper("clipRange2", "-CT 1-5,11-15", "be4fcad5b666a5540028b774169cbad7", "ec4cf54ed50a6baf69dbf98782c19aeb"); }

    @Test public void testClipSeq() { testClipper("clipSeqX", "-X CCCCC", "db199bd06561c9f2122f6ffb07941fbc", "a9cf540e4ed2514061248a878e09a09c"); }
    @Test public void testClipSeqFile() { testClipper("clipSeqXF", "-XF " + privateTestDir + "seqsToClip.fasta", "d011a3152b31822475afbe0281491f8d", "906871df304dd966682e5798d59fc86b"); }

    @Test public void testClipMulti() { testClipper("clipSeqMulti", "-QT 10 -CT 1-5 -XF " + privateTestDir + "seqsToClip.fasta -X CCCCC", "a23187bd9bfb06557f799706d98441de", "b41995fea04034ca0427c4a71504ef83"); }

    @Test public void testClipNs() { testClipper("testClipNs", "-QT 10 -CR WRITE_NS", Q10ClipOutput, "27847d330b962e60650df23b6efc8c3c"); }
    @Test public void testClipQ0s() { testClipper("testClipQs", "-QT 10 -CR WRITE_Q0S", Q10ClipOutput, "195b8bdfc0186fdca742764aa9b06363"); }
    @Test public void testClipSoft() { testClipper("testClipSoft", "-QT 10 -CR SOFTCLIP_BASES", Q10ClipOutput, "08d16051be0b3fa3453eb1e6ca48b098"); }

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
                Arrays.asList("a2819d54b2110150e38511f5a55db91d", "55c01ccc2e84481b22d3632cdb06c8ba"));
        executeTest("clipOriginalQuals", spec);
    }
}
