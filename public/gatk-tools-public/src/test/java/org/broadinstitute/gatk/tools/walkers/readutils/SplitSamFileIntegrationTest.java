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
import java.util.Collections;

public class SplitSamFileIntegrationTest extends WalkerTest {

    @Test
    public void testSplitSamFile() {
        final String prefix = "splitsam";
        WalkerTestSpec spec = new WalkerTestSpec(
                "-T SplitSamFile" +
                        " -R " + b37KGReference +
                        " -I " + privateTestDir+"/CEUTrio.HiSeq.b37.MT.1_50.bam" +
                        " --outputRoot " + prefix,
                Collections.<String>emptyList()
        );
        addSplitOutput(spec, prefix, "NA12878", "356bddf614a44a0234f7b116f5cde998");
        addSplitOutput(spec, prefix, "NA12891", "5dd92a0089e3f781e77d94010264b9fd");
        addSplitOutput(spec, prefix, "NA12892", "ea8f16c18c89b49ce9d32d76c6af75f9");
        executeTest("testSplitSamFile", spec);
    }

    private void addSplitOutput(final WalkerTestSpec spec, final String outputPrefix, final String sample, final String md5) {
        final File outputFile = new File(outputPrefix + sample + ".bam");
        spec.addAuxFile(md5, outputFile);

        //The AuxFile mechanism will ensure the bam is deleted, but it doesn't know about indices
        new File(outputFile.getAbsolutePath() + ".bai").deleteOnExit();
        new File(outputFile.getAbsolutePath().replaceAll("bam$", ".bai")).deleteOnExit();
    }
}
