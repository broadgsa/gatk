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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Test basic functionality in SeekableBufferedStream.
 */
public class SeekableBufferedStreamUnitTest extends BaseTest {
    private static File InputFile = new File(validationDataLocation + "megabyteZeros.dat");

    final private int BUFFERED_STREAM_BUFFER_SIZE = 100;
    private byte buffer[] = new byte[BUFFERED_STREAM_BUFFER_SIZE * 10];


    @DataProvider(name = "BasicArgumentsDivisible")
    public Integer[][] DivisableReads() {
        return new Integer[][]{{1}, {4}, {5}, {10}, {20}, {50}, {100}};
    }

    @DataProvider(name = "BasicArgumentsIndivisibleAndSmall")
    public Integer[][] InDivisableReadsSmall() {
        return new Integer[][]{{3}, {11}, {31}, {51}, {77}, {99}};
    }

    @DataProvider(name = "BasicArgumentsIndivisibleYetLarge")
    public Integer[][] InDivisableReadsLarge() {
        return new Integer[][]{{101}, {151}, {205}, {251}, {301}};
    }


    private void testReadsLength(int length) throws IOException {
        final int READ_SIZE=100000; //file is 10^6, so make this smaller to be safe.

        SeekableFileStream fileStream = new SeekableFileStream(InputFile);
        SeekableBufferedStream bufferedStream = new SeekableBufferedStream(fileStream, BUFFERED_STREAM_BUFFER_SIZE);

        for (int i = 0; i < READ_SIZE / length; ++i) {
            Assert.assertEquals(bufferedStream.read(buffer, 0, length), length);
        }

    }

    // These tests fail because SeekableBuffered stream may return _less_ than the amount you are asking for.
    // make sure that you wrap reads with while-loops.  If these test start failing (meaning that the reads work properly,
    // the layer of protection built into GATKBamIndex can be removed.
    //
    // pdexheimer, Jan 2015 - SeekableBufferedStream no longer returns less than the expected amount.
    // Renaming testIndivisableSmallReadsFAIL to testIndivisableSmallReadsPASS and removing the expected exception
    // If this bug regresses, the while loop will need to be re-introduced into GATKBamIndex.read()

    @Test(dataProvider = "BasicArgumentsIndivisibleAndSmall", enabled = true)
    public void testIndivisableSmallReadsPASS(Integer readLength) throws IOException {
        testReadsLength(readLength);
    }

    //Evidently, if you ask for a read length that's larger than the inernal buffer,
    //SeekableBufferedStreamdoes something else and gives you what you asked for

    @Test(dataProvider = "BasicArgumentsIndivisibleYetLarge", enabled = true)
    public void testIndivisableLargeReadsPASS(Integer readLength) throws IOException {
        testReadsLength(readLength);
    }

    // if the readlength divides the buffer, there are no failures
    @Test(dataProvider = "BasicArgumentsDivisible", enabled = true)
    public void testDivisableReadsPASS(Integer readLength) throws IOException {
        testReadsLength(readLength);
    }


}
