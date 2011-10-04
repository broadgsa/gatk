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
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;

/**
 * Test the GATK core interval parsing mechanism.
 */
public class IntervalIntegrationTest extends WalkerTest {
    @Test(enabled = true)
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

    @Test(enabled = true)
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

    @Test
    public void testUnmappedReadInclusion() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T PrintReads" +
                        " -I " + validationDataLocation + "MV1994.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -L unmapped" +
                        " -U",
                        0, // two output files
                        Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadInclusion",".bam");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("748a38ed5eb0a043dfc7b82f0d1e8063",createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("fadcdf88597b9609c5f2a17f4c6eb455", createTempFileFromBase(baseOutputFile.getAbsolutePath().substring(0,baseOutputFile.getAbsolutePath().indexOf(".bam"))+".bai"));

        executeTest("testUnmappedReadInclusion",spec);
    }

    @Test(enabled = true)
    public void testUnmappedReadExclusion() {
        WalkerTest.WalkerTestSpec spec = new WalkerTest.WalkerTestSpec(
                "-T PrintReads" +
                        " -I " + validationDataLocation + "MV1994.bam" +
                        " -R " + validationDataLocation + "Escherichia_coli_K12_MG1655.fasta" +
                        " -XL unmapped" +
                        " -U",
                        0, // two output files
                        Collections.<String>emptyList());

        // our base file
        File baseOutputFile = createTempFile("testUnmappedReadExclusion",".bam");
        spec.setOutputFileLocation(baseOutputFile);
        spec.addAuxFile("8236f0b2df5a692e54751b08bc3836fa",createTempFileFromBase(baseOutputFile.getAbsolutePath()));
        spec.addAuxFile("b341d808ecc33217f37c0c0cde2a3e2f", createTempFileFromBase(baseOutputFile.getAbsolutePath().substring(0,baseOutputFile.getAbsolutePath().indexOf(".bam"))+".bai"));

        executeTest("testUnmappedReadExclusion",spec);
    }


}
