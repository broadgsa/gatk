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

package org.broadinstitute.gatk.engine.datasources.reference;

import org.broadinstitute.gatk.engine.walkers.WalkerTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;

public class ReferenceDataSourceIntegrationTest extends WalkerTest {

    @Test
    public void testReferenceWithMissingFaiFile() throws IOException {
        final File dummyReference = createTempFile("dummy", ".fasta");
        final File dictFile = new File(dummyReference.getAbsolutePath().replace(".fasta", ".dict"));
        dictFile.deleteOnExit();
        Assert.assertTrue(dictFile.createNewFile());

        final WalkerTestSpec spec = new WalkerTestSpec(
            " -T TestPrintReadsWalker" +
            " -R " + dummyReference.getAbsolutePath() +
            " -I " + privateTestDir + "NA12878.4.snippet.bam" +
            " -o %s",
            1,
            UserException.MissingReferenceFaiFile.class
        );

        executeTest("testReferenceWithMissingFaiFile", spec);
    }

    @Test
    public void testReferenceWithMissingDictFile() throws IOException {
        final File dummyReference = createTempFile("dummy", ".fasta");
        final File faiFile = new File(dummyReference.getAbsolutePath() + ".fai");
        faiFile.deleteOnExit();
        Assert.assertTrue(faiFile.createNewFile());

        final WalkerTestSpec spec = new WalkerTestSpec(
                " -T TestPrintReadsWalker" +
                " -R " + dummyReference.getAbsolutePath() +
                " -I " + privateTestDir + "NA12878.4.snippet.bam" +
                " -o %s",
                1,
                UserException.MissingReferenceDictFile.class
        );

        executeTest("testReferenceWithMissingDictFile", spec);
    }
}
