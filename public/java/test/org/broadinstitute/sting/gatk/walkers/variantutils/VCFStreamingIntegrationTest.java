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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.WalkerTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;

/**
 * Test the generic VCF streaming workflow with a few basic VCF routines.
 */
public class VCFStreamingIntegrationTest extends WalkerTest {
    @Test
    public void testSimpleVCFStreaming() throws IOException {
        // Create a FIFO.  This seems to be the only way to create an interprocess FIFO in Java (java.nio.Pipe is intraprocess only).
        File tmpFifo = File.createTempFile("vcfstreaming","");
        Runtime.getRuntime().exec(new String[] {"mkfifo",tmpFifo.getAbsolutePath()});


        // Copy VCF data from the test file into the FIFO.
        String testFile = validationDataLocation + "yri.trio.gatk.ug.head.vcf";
        FileInputStream inputStream = new FileInputStream(testFile);
        FileOutputStream outputStream = new FileOutputStream(tmpFifo);
        outputStream.getChannel().transferFrom(inputStream.getChannel(),0,inputStream.getChannel().size());
        outputStream.close();
        inputStream.close();

        WalkerTestSpec spec = new WalkerTestSpec(
            "-T SelectVariants" +
                    " -R " + b36KGReference +
                    " --variant:vcf3,storage=STREAM " + tmpFifo.getAbsolutePath() +
                    " --NO_HEADER" +
                    " -o %s",
            1,
            Arrays.asList("658f580f7a294fd334bd897102616fed")
        );

        executeTest("testSimpleVCFStreaming", spec);

        tmpFifo.delete();
    }

    @Test
    public void testVCFStreamingChain() throws IOException {
        // Create a FIFO.  This seems to be the only way to create an interprocess FIFO in Java (java.nio.Pipe is intraprocess only).
        File tmpFifo = File.createTempFile("vcfstreaming","");
        Runtime.getRuntime().exec(new String[] {"mkfifo",tmpFifo.getAbsolutePath()});

        String testFile = validationDataLocation + "yri.trio.gatk.ug.head.vcf";

        // Output select to FIFO
        WalkerTestSpec selectTestSpec = new WalkerTestSpec(
            "-T SelectVariants" +
            " -R " + b36KGReference +
            " --variant:vcf3,storage=STREAM " + testFile +
            " --NO_HEADER" +
            " -select 'QD > 2.0'" +
            " -o " + tmpFifo.getAbsolutePath(),
            0,
            Collections.<String>emptyList()
        );
        executeTest("testVCFStreamingChain", selectTestSpec);

        // Eval compare the full set to the subselection.
        selectTestSpec = new WalkerTestSpec(
            "-T VariantEval" +
            " -R " + b36KGReference +
            " --eval:vcf3 " + testFile +
            " --comp:vcf,storage=STREAM " + tmpFifo.getAbsolutePath() +
            " -EV CompOverlap -noEV -noST" +
            " -o %s",
            1,
            Arrays.asList("d46a735ffa898f4aa6b3758c5b03f06d")
        );
        executeTest("testVCFStreamingChain", selectTestSpec);

        tmpFifo.delete();
    }
}
