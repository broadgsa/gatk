/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.downsampling;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSingleSampleReadStream;
import org.broadinstitute.sting.utils.sam.ArtificialSingleSampleReadStreamAnalyzer;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public class DownsamplingReadsIteratorUnitTest extends BaseTest {

    private static class DownsamplingReadsIteratorTest extends TestDataProvider {
        private DownsamplingReadsIterator downsamplingIter;
        private int targetCoverage;
        private ArtificialSingleSampleReadStream stream;
        private ArtificialSingleSampleReadStreamAnalyzer streamAnalyzer;

        public DownsamplingReadsIteratorTest( ArtificialSingleSampleReadStream stream, int targetCoverage ) {
            super(DownsamplingReadsIteratorTest.class);

            this.stream = stream;
            this.targetCoverage = targetCoverage;

            setName(String.format("%s: targetCoverage=%d numContigs=%d stacksPerContig=%d readsPerStack=%d-%d distanceBetweenStacks=%d-%d readLength=%d-%d unmappedReads=%d",
                    getClass().getSimpleName(),
                    targetCoverage,
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
            streamAnalyzer = new PositionallyDownsampledArtificialSingleSampleReadStreamAnalyzer(stream, targetCoverage);
            downsamplingIter = new DownsamplingReadsIterator(stream.getStingSAMIterator(), new SimplePositionalDownsampler<SAMRecord>(targetCoverage));

            streamAnalyzer.analyze(downsamplingIter);

            // Check whether the observed properties of the downsampled stream are what they should be
            streamAnalyzer.validate();

            // Allow memory used by this test to be reclaimed
            stream = null;
            streamAnalyzer = null;
            downsamplingIter = null;
        }
    }

    @DataProvider(name = "DownsamplingReadsIteratorTestDataProvider")
    public Object[][] createDownsamplingReadsIteratorTests() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(5, 1, 10000);
        String readGroupID = "testReadGroup";
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupID);
        readGroup.setSample("testSample");
        header.addReadGroup(readGroup);

        // Values that don't vary across tests
        int targetCoverage = 10;
        int minReadLength = 50;
        int maxReadLength = 100;
        int minDistanceBetweenStacks = 1;
        int maxDistanceBetweenStacks = maxReadLength + 1;

        GenomeAnalysisEngine.resetRandomGenerator();

        // brute force testing!
        for ( int numContigs : Arrays.asList(1, 2, 5) ) {
            for ( int stacksPerContig : Arrays.asList(1, 2, 10) ) {
                for ( int minReadsPerStack : Arrays.asList(1, targetCoverage / 2, targetCoverage, targetCoverage - 1, targetCoverage + 1, targetCoverage * 2) ) {
                    for ( int maxReadsPerStack : Arrays.asList(1, targetCoverage / 2, targetCoverage, targetCoverage - 1, targetCoverage + 1, targetCoverage * 2) ) {
                        for ( int numUnmappedReads : Arrays.asList(0, 1, targetCoverage, targetCoverage * 2) ) {
                            // Only interested in sane read stream configurations here
                            if ( minReadsPerStack <= maxReadsPerStack ) {
                                new DownsamplingReadsIteratorTest(new ArtificialSingleSampleReadStream(header,
                                                                                                       readGroupID,
                                                                                                       numContigs,
                                                                                                       stacksPerContig,
                                                                                                       minReadsPerStack,
                                                                                                       maxReadsPerStack,
                                                                                                       minDistanceBetweenStacks,
                                                                                                       maxDistanceBetweenStacks,
                                                                                                       minReadLength,
                                                                                                       maxReadLength,
                                                                                                       numUnmappedReads),
                                                                  targetCoverage);
                            }
                        }
                    }
                }
            }
        }

        return DownsamplingReadsIteratorTest.getTests(DownsamplingReadsIteratorTest.class);
    }

    @Test(dataProvider = "DownsamplingReadsIteratorTestDataProvider")
    public void runDownsamplingReadsIteratorTest( DownsamplingReadsIteratorTest test ) {
        logger.warn("Running test: " + test);

        GenomeAnalysisEngine.resetRandomGenerator();
        test.run();
    }
}
