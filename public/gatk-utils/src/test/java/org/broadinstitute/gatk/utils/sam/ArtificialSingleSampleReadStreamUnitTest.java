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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import org.broadinstitute.gatk.utils.BaseTest;

public class ArtificialSingleSampleReadStreamUnitTest extends BaseTest {

    private static class ArtificialSingleSampleReadStreamTest extends TestDataProvider {
        private ArtificialSingleSampleReadStream stream;
        private ArtificialSingleSampleReadStreamAnalyzer streamAnalyzer;

        public ArtificialSingleSampleReadStreamTest( ArtificialSingleSampleReadStream stream ) {
            super(ArtificialSingleSampleReadStreamTest.class);

            this.stream = stream;

            setName(String.format("%s: numContigs=%d stacksPerContig=%d readsPerStack=%d-%d distanceBetweenStacks=%d-%d readLength=%d-%d unmappedReads=%d",
                    getClass().getSimpleName(),
                    stream.getNumContigs(),
                    stream.getNumStacksPerContig(),
                    stream.getMinReadsPerStack(),
                    stream.getMaxReadsPerStack(),
                    stream.getMinDistanceBetweenStacks(),
                    stream.getMaxDistanceBetweenStacks(),
                    stream.getMinReadLength(),
                    stream.getMaxReadLength(),
                    stream.getNumUnmappedReads()));
        }

        public void run() {
            streamAnalyzer= new ArtificialSingleSampleReadStreamAnalyzer(stream);

            streamAnalyzer.analyze(stream);

            // Check whether the observed properties of the stream match its nominal properties
            streamAnalyzer.validate();
        }
    }

    @DataProvider(name = "ArtificialSingleSampleReadStreamTestDataProvider")
    public Object[][] createArtificialSingleSampleReadStreamTests() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, 10000);
        String readGroupID = "testReadGroup";
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupID);
        readGroup.setSample("testSample");
        header.addReadGroup(readGroup);

        Utils.resetRandomGenerator();

        // brute force testing!
        for ( int numContigs = 0; numContigs <= 2; numContigs++ ) {
            for ( int stacksPerContig = 0; stacksPerContig <= 2; stacksPerContig++ ) {
                for ( int minReadsPerStack = 1; minReadsPerStack <= 2; minReadsPerStack++ ) {
                    for ( int maxReadsPerStack = 1; maxReadsPerStack <= 3; maxReadsPerStack++ ) {
                        for ( int minDistanceBetweenStacks = 1; minDistanceBetweenStacks <= 2; minDistanceBetweenStacks++ ) {
                            for ( int maxDistanceBetweenStacks = 1; maxDistanceBetweenStacks <= 3; maxDistanceBetweenStacks++ ) {
                                for ( int minReadLength = 1; minReadLength <= 2; minReadLength++ ) {
                                    for ( int maxReadLength = 1; maxReadLength <= 3; maxReadLength++ ) {
                                        for ( int numUnmappedReads = 0; numUnmappedReads <= 2; numUnmappedReads++ ) {
                                            // Only test sane combinations here
                                            if ( minReadsPerStack <= maxReadsPerStack &&
                                                 minDistanceBetweenStacks <= maxDistanceBetweenStacks &&
                                                 minReadLength <= maxReadLength &&
                                                 ((numContigs > 0 && stacksPerContig > 0) || (numContigs == 0 && stacksPerContig == 0)) ) {

                                                new ArtificialSingleSampleReadStreamTest(new ArtificialSingleSampleReadStream(header,
                                                                                                                              readGroupID,
                                                                                                                              numContigs,
                                                                                                                              stacksPerContig,
                                                                                                                              minReadsPerStack,
                                                                                                                              maxReadsPerStack,
                                                                                                                              minDistanceBetweenStacks,
                                                                                                                              maxDistanceBetweenStacks,
                                                                                                                              minReadLength,
                                                                                                                              maxReadLength,
                                                                                                                              numUnmappedReads));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return ArtificialSingleSampleReadStreamTest.getTests(ArtificialSingleSampleReadStreamTest.class);
    }

    @Test(dataProvider = "ArtificialSingleSampleReadStreamTestDataProvider")
    public void testArtificialSingleSampleReadStream( ArtificialSingleSampleReadStreamTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();
        test.run();
    }

    @DataProvider(name = "ArtificialSingleSampleReadStreamInvalidArgumentsTestDataProvider")
    public Object[][] createInvalidArgumentsTests() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(3, 1, 10000);
        String readGroupID = "testReadGroup";
        header.addReadGroup(new SAMReadGroupRecord(readGroupID));

        return new Object[][] {
            {"testNullHeader", null, readGroupID, 1, 1, 1, 2, 1, 2, 1, 2, 0},
            {"testNullReadGroup", header, null, 1, 1, 1, 2, 1, 2, 1, 2, 0},
            {"testInvalidReadGroup", header, "foo", 1, 1, 1, 2, 1, 2, 1, 2, 0},
            {"testInvalidNumContigs", header, readGroupID, -1, 1, 1, 2, 1, 2, 1, 2, 0},
            {"testInvalidNumStacksPerContig", header, readGroupID, 1, -1, 1, 2, 1, 2, 1, 2, 0},
            {"test0ContigsNon0StacksPerContig", header, readGroupID, 0, 1, 1, 2, 1, 2, 1, 2, 0},
            {"testNon0Contigs0StacksPerContig", header, readGroupID, 1, 0, 1, 2, 1, 2, 1, 2, 0},
            {"testInvalidMinReadsPerStack", header, readGroupID, 1, 1, -1, 2, 1, 2, 1, 2, 0},
            {"testInvalidMaxReadsPerStack", header, readGroupID, 1, 1, 1, -2, 1, 2, 1, 2, 0},
            {"testInvalidMinDistanceBetweenStacks", header, readGroupID, 1, 1, 1, 2, -1, 2, 1, 2, 0},
            {"testInvalidMaxDistanceBetweenStacks", header, readGroupID, 1, 1, 1, 2, 1, -2, 1, 2, 0},
            {"testInvalidMinReadLength", header, readGroupID, 1, 1, 1, 2, 1, 2, -1, 2, 0},
            {"testInvalidMaxReadLength", header, readGroupID, 1, 1, 1, 2, 1, 2, 1, -2, 0},
            {"testInvalidReadsPerStackRange", header, readGroupID, 1, 1, 2, 1, 1, 2, 1, 2, 0},
            {"testInvalidDistanceBetweenStacksRange", header, readGroupID, 1, 1, 1, 2, 2, 1, 1, 2, 0},
            {"testInvalidReadLengthRange", header, readGroupID, 1, 1, 1, 2, 1, 2, 2, 1, 0},
            {"testInvalidNumUnmappedReads", header, readGroupID, 1, 1, 1, 2, 1, 2, 1, 2, -1},
        };
    }

    @Test(dataProvider = "ArtificialSingleSampleReadStreamInvalidArgumentsTestDataProvider",
          expectedExceptions = ReviewedGATKException.class)
    public void testInvalidArguments( String testName,
                                      SAMFileHeader header,
                                      String readGroupID,
                                      int numContigs,
                                      int numStacksPerContig,
                                      int minReadsPerStack,
                                      int maxReadsPerStack,
                                      int minDistanceBetweenStacks,
                                      int maxDistanceBetweenStacks,
                                      int minReadLength,
                                      int maxReadLength,
                                      int numUnmappedReads ) {

        logger.warn("Running test: " + testName);

        ArtificialSingleSampleReadStream stream = new ArtificialSingleSampleReadStream(header,
                                                                                       readGroupID,
                                                                                       numContigs,
                                                                                       numStacksPerContig,
                                                                                       minReadsPerStack,
                                                                                       maxReadsPerStack,
                                                                                       minDistanceBetweenStacks,
                                                                                       maxDistanceBetweenStacks,
                                                                                       minReadLength,
                                                                                       maxReadLength,
                                                                                       numUnmappedReads);
    }
}
