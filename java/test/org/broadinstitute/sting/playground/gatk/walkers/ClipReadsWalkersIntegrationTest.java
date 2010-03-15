package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.WalkerTest;
import org.junit.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ClipReadsWalkersIntegrationTest extends WalkerTest {
    public void testClipper(String name, String args, String md51, String md52) {
        WalkerTestSpec spec = new WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta " +
                        "-T ClipReads " +
                        "-I " + validationDataLocation + "clippingReadsTest.bam " +
                        "-o %s " +
                        "-ob %s " + args,
                2, // just one output file
                Arrays.asList("tmp", "bam"),
                Arrays.asList(md51, md52));
        List<File> result = executeTest(name, spec).getFirst();
    }

    final static String Q10ClipOutput = "b29c5bc1cb9006ed9306d826a11d444f";
    @Test public void testQClip0() { testClipper("clipQSum0", "-QT 0", "117a4760b54308f81789c39b1c9de578", "2465660bcd975a1dc6dfbf40a21bf6ad"); }
    @Test public void testQClip2() { testClipper("clipQSum2", "-QT 2", Q10ClipOutput, "fb77d3122df468a71e03ca92b69493f4"); }
    @Test public void testQClip10() { testClipper("clipQSum10", "-QT 10", "b29c5bc1cb9006ed9306d826a11d444f", "fb77d3122df468a71e03ca92b69493f4"); }
    @Test public void testQClip20() { testClipper("clipQSum20", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "9a4b1c83c026ca83db00bb71999246cf"); }
    @Test public void testQClip30() { testClipper("clipQSum30", "-QT 20", "6c3434dce66ae5c9eeea502f10fb9bee", "9a4b1c83c026ca83db00bb71999246cf"); }

    @Test public void testClipRange1() { testClipper("clipRange1", "-CT 1-5", "b5acd753226e25b1e088838c1aab9117", "9a08474a13fbb897b7c9dca58d19884f"); }
    @Test public void testClipRange2() { testClipper("clipRange2", "-CT 1-5,11-15", "be4fcad5b666a5540028b774169cbad7", "f05ab5fe821b77cd5b066212ff56f8ff"); }

    @Test public void testClipSeq() { testClipper("clipSeqX", "-X CCCCC", "db199bd06561c9f2122f6ffb07941fbc", "c218c0649838423a06f3296430f65c4f"); }
    @Test public void testClipSeqFile() { testClipper("clipSeqXF", "-XF " + validationDataLocation + "seqsToClip.fasta", "d011a3152b31822475afbe0281491f8d", "1151e10833da794203df2ba7cc76d5c5"); }

    @Test public void testClipMulti() { testClipper("clipSeqMulti", "-QT 10 -CT 1-5 -XF " + validationDataLocation + "seqsToClip.fasta -X CCCCC", "a23187bd9bfb06557f799706d98441de", "4a1153d6f0600cf53ff7959a043e57cc"); }

    @Test public void testClipNs() { testClipper("testClipNs", "-QT 10 -CR WRITE_NS", Q10ClipOutput, "fb77d3122df468a71e03ca92b69493f4"); }
    @Test public void testClipQ0s() { testClipper("testClipQs", "-QT 10 -CR WRITE_Q0S", Q10ClipOutput, "24053a87b00c0bc2ddf420975e9fea4d"); }
    @Test public void testClipSoft() { testClipper("testClipSoft", "-QT 10 -CR SOFTCLIP_BASES", Q10ClipOutput, "aeb67cca75285a68af8a965faa547e7f"); }

    @Test
    public void testUseOriginalQuals() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-R " + seqLocation + "references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta" +
                        " -T ClipReads" +
                        " -I " + validationDataLocation + "originalQuals.chr1.1-1K.bam" +
                        " -L chr1:1-1,000" +
                        " -OQ -QT 4" +
                        " -o %s -ob %s",
                2,
                Arrays.asList("55c01ccc2e84481b22d3632cdb06c8ba", "f9b1347fabbc33bb24f7c7fa8dfb798b"));
        executeTest("clipOriginalQuals", spec);
    }
}