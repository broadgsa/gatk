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

package org.broadinstitute.gatk.engine.downsampling;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.downsampling.PerSampleDownsamplingReadsIterator;
import org.broadinstitute.gatk.utils.downsampling.ReadsDownsamplerFactory;
import org.broadinstitute.gatk.utils.downsampling.SimplePositionalDownsamplerFactory;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.engine.iterators.VerifyingSamIterator;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.sam.ArtificialMultiSampleReadStream;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.ArtificialSingleSampleReadStream;
import org.broadinstitute.gatk.utils.sam.ArtificialSingleSampleReadStreamAnalyzer;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class PerSampleDownsamplingReadsIteratorUnitTest extends BaseTest {

    private static class PerSampleDownsamplingReadsIteratorTest extends TestDataProvider {

        // TODO: tests should distinguish between variance across samples and variance within a sample

        private enum StreamDensity {
            SPARSE         (MAX_READ_LENGTH,     MAX_READ_LENGTH * 2),
            DENSE          (1,                   MIN_READ_LENGTH),
            MIXED          (1,                   MAX_READ_LENGTH * 2),
            UNIFORM_DENSE  (1,                   1),
            UNIFORM_SPARSE (MAX_READ_LENGTH * 2, MAX_READ_LENGTH * 2);

            int minDistanceBetweenStacks;
            int maxDistanceBetweenStacks;

            StreamDensity( int minDistanceBetweenStacks, int maxDistanceBetweenStacks ) {
                this.minDistanceBetweenStacks = minDistanceBetweenStacks;
                this.maxDistanceBetweenStacks = maxDistanceBetweenStacks;
            }

            public String toString() {
                return String.format("StreamDensity:%d-%d", minDistanceBetweenStacks, maxDistanceBetweenStacks);
            }
        }

        private enum StreamStackDepth {
            NON_UNIFORM_LOW   (1,  5),
            NON_UNIFORM_HIGH  (15, 20),
            NON_UNIFORM_MIXED (1,  20),
            UNIFORM_SINGLE    (1,  1),
            UNIFORM_LOW       (2,  2),
            UNIFORM_HIGH      (20, 20),
            UNIFORM_MEDIUM    (10, 10);   // should set target coverage to this value for testing

            int minReadsPerStack;
            int maxReadsPerStack;

            StreamStackDepth( int minReadsPerStack, int maxReadsPerStack ) {
                this.minReadsPerStack = minReadsPerStack;
                this.maxReadsPerStack = maxReadsPerStack;
            }

            public boolean isUniform() {
                return minReadsPerStack == maxReadsPerStack;
            }

            public String toString() {
                return String.format("StreamStackDepth:%d-%d", minReadsPerStack, maxReadsPerStack);
            }
        }

        private enum StreamStacksPerContig {
            UNIFORM(20, 20),
            NON_UNIFORM(1, 30);

            int minStacksPerContig;
            int maxStacksPerContig;

            StreamStacksPerContig( int minStacksPerContig, int maxStacksPerContig ) {
                this.minStacksPerContig = minStacksPerContig;
                this.maxStacksPerContig = maxStacksPerContig;
            }

            public boolean isUniform() {
                return minStacksPerContig == maxStacksPerContig;
            }

            public String toString() {
                return String.format("StreamStacksPerContig:%d-%d", minStacksPerContig, maxStacksPerContig);
            }
        }

        // Not interested in testing multiple ranges for the read lengths, as none of our current
        // downsamplers are affected by read length
        private static final int MIN_READ_LENGTH = 50;
        private static final int MAX_READ_LENGTH = 150;

        private ReadsDownsamplerFactory<SAMRecord> downsamplerFactory;
        private int targetCoverage;
        private int numSamples;
        private int minContigs;
        private int maxContigs;
        private StreamDensity streamDensity;
        private StreamStackDepth streamStackDepth;
        private StreamStacksPerContig streamStacksPerContig;
        private double unmappedReadsFraction;
        private int unmappedReadsCount;
        private boolean verifySortedness;

        private ArtificialMultiSampleReadStream mergedReadStream;
        private Map<String, ArtificialSingleSampleReadStream> perSampleArtificialReadStreams;
        private Map<String, ArtificialSingleSampleReadStreamAnalyzer> perSampleStreamAnalyzers;
        private SAMFileHeader header;

        public PerSampleDownsamplingReadsIteratorTest( ReadsDownsamplerFactory<SAMRecord> downsamplerFactory,
                                                       int targetCoverage,
                                                       int numSamples,
                                                       int minContigs,
                                                       int maxContigs,
                                                       StreamDensity streamDensity,
                                                       StreamStackDepth streamStackDepth,
                                                       StreamStacksPerContig streamStacksPerContig,
                                                       double unmappedReadsFraction,
                                                       int unmappedReadsCount,
                                                       boolean verifySortedness ) {
            super(PerSampleDownsamplingReadsIteratorTest.class);

            this.downsamplerFactory = downsamplerFactory;
            this.targetCoverage = targetCoverage;
            this.numSamples = numSamples;
            this.minContigs = minContigs;
            this.maxContigs = maxContigs;
            this.streamDensity = streamDensity;
            this.streamStackDepth = streamStackDepth;
            this.streamStacksPerContig = streamStacksPerContig;
            this.unmappedReadsFraction = unmappedReadsFraction;
            this.unmappedReadsCount = unmappedReadsCount;
            this.verifySortedness = verifySortedness;

            header = createHeader();
            createReadStreams();

            setName(String.format("%s: targetCoverage=%d numSamples=%d minContigs=%d maxContigs=%d %s %s %s unmappedReadsFraction=%.2f unmappedReadsCount=%d verifySortedness=%b",
                    getClass().getSimpleName(), targetCoverage, numSamples, minContigs, maxContigs, streamDensity, streamStackDepth, streamStacksPerContig, unmappedReadsFraction, unmappedReadsCount, verifySortedness));
        }

        private SAMFileHeader createHeader() {
            SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(maxContigs, 1, (streamDensity.maxDistanceBetweenStacks + MAX_READ_LENGTH) * streamStacksPerContig.maxStacksPerContig + 100000);
            List<String> readGroups = new ArrayList<String>(numSamples);
            List<String> sampleNames = new ArrayList<String>(numSamples);

            for ( int i = 0; i < numSamples; i++ ) {
                readGroups.add("ReadGroup" + i);
                sampleNames.add("Sample" + i);
            }

            return ArtificialSAMUtils.createEnumeratedReadGroups(header, readGroups, sampleNames);
        }

        private void createReadStreams() {
            perSampleArtificialReadStreams = new HashMap<String, ArtificialSingleSampleReadStream>(numSamples);
            perSampleStreamAnalyzers = new HashMap<String, ArtificialSingleSampleReadStreamAnalyzer>(numSamples);

            for (SAMReadGroupRecord readGroup : header.getReadGroups() ) {
                String readGroupID = readGroup.getReadGroupId();
                String sampleName = readGroup.getSample();

                int thisSampleNumContigs = MathUtils.randomIntegerInRange(minContigs, maxContigs);
                int thisSampleStacksPerContig = MathUtils.randomIntegerInRange(streamStacksPerContig.minStacksPerContig, streamStacksPerContig.maxStacksPerContig);

                int thisSampleNumUnmappedReads = Utils.getRandomGenerator().nextDouble() < unmappedReadsFraction ? unmappedReadsCount : 0;

                ArtificialSingleSampleReadStream thisSampleStream = new ArtificialSingleSampleReadStream(header,
                                                                                                         readGroupID,
                                                                                                         thisSampleNumContigs,
                                                                                                         thisSampleStacksPerContig,
                                                                                                         streamStackDepth.minReadsPerStack,
                                                                                                         streamStackDepth.maxReadsPerStack,
                                                                                                         streamDensity.minDistanceBetweenStacks,
                                                                                                         streamDensity.maxDistanceBetweenStacks,
                                                                                                         MIN_READ_LENGTH,
                                                                                                         MAX_READ_LENGTH,
                                                                                                         thisSampleNumUnmappedReads);
                perSampleArtificialReadStreams.put(sampleName, thisSampleStream);
                perSampleStreamAnalyzers.put(sampleName, new PositionallyDownsampledArtificialSingleSampleReadStreamAnalyzer(thisSampleStream, targetCoverage));
            }

            mergedReadStream = new ArtificialMultiSampleReadStream(perSampleArtificialReadStreams.values());
        }

        public void run() {
            GATKSAMIterator downsamplingIter = new PerSampleDownsamplingReadsIterator(mergedReadStream.getGATKSAMIterator(), downsamplerFactory);

            if ( verifySortedness ) {
                downsamplingIter = new VerifyingSamIterator(downsamplingIter);
            }

            while ( downsamplingIter.hasNext() ) {
                SAMRecord read = downsamplingIter.next();
                String sampleName = read.getReadGroup() != null ? read.getReadGroup().getSample() : null;

                ArtificialSingleSampleReadStreamAnalyzer analyzer = perSampleStreamAnalyzers.get(sampleName);
                if ( analyzer != null ) {
                    analyzer.update(read);
                }
                else {
                    throw new ReviewedGATKException("bug: stream analyzer for sample " + sampleName + " not found");
                }
            }

            for ( Map.Entry<String, ArtificialSingleSampleReadStreamAnalyzer> analyzerEntry : perSampleStreamAnalyzers.entrySet() ) {
                ArtificialSingleSampleReadStreamAnalyzer analyzer = analyzerEntry.getValue();
                analyzer.finalizeStats();

                // Validate the downsampled read stream for each sample individually
                analyzer.validate();
            }

            // Allow memory used by this test to be reclaimed:
            mergedReadStream = null;
            perSampleArtificialReadStreams = null;
            perSampleStreamAnalyzers = null;
        }
    }

    @DataProvider(name = "PerSampleDownsamplingReadsIteratorTestDataProvider")
    public Object[][] createPerSampleDownsamplingReadsIteratorTests() {

        Utils.resetRandomGenerator();

        // Some values don't vary across tests
        int targetCoverage = PerSampleDownsamplingReadsIteratorTest.StreamStackDepth.UNIFORM_MEDIUM.minReadsPerStack;
        ReadsDownsamplerFactory<SAMRecord> downsamplerFactory = new SimplePositionalDownsamplerFactory<SAMRecord>(targetCoverage);
        int maxContigs = 3;
        boolean verifySortedness = true;

        for ( int numSamples : Arrays.asList(1, 2, 10) ) {
            for ( int minContigs = 1; minContigs <= maxContigs; minContigs++ ) {
                for ( PerSampleDownsamplingReadsIteratorTest.StreamDensity streamDensity : PerSampleDownsamplingReadsIteratorTest.StreamDensity.values() ) {
                    for ( PerSampleDownsamplingReadsIteratorTest.StreamStackDepth streamStackDepth : PerSampleDownsamplingReadsIteratorTest.StreamStackDepth.values() ) {
                        for (PerSampleDownsamplingReadsIteratorTest.StreamStacksPerContig streamStacksPerContig : PerSampleDownsamplingReadsIteratorTest.StreamStacksPerContig.values() ) {
                            for ( double unmappedReadsFraction : Arrays.asList(0.0, 1.0, 0.5) ) {
                                for ( int unmappedReadsCount : Arrays.asList(1, 50) ) {
                                    new PerSampleDownsamplingReadsIteratorTest(downsamplerFactory,
                                                                               targetCoverage,
                                                                               numSamples,
                                                                               minContigs,
                                                                               maxContigs,
                                                                               streamDensity,
                                                                               streamStackDepth,
                                                                               streamStacksPerContig,
                                                                               unmappedReadsFraction,
                                                                               unmappedReadsCount,
                                                                               verifySortedness);
                                }
                            }
                        }
                    }
                }
            }
        }

        return PerSampleDownsamplingReadsIteratorTest.getTests(PerSampleDownsamplingReadsIteratorTest.class);
    }

    @Test(dataProvider = "PerSampleDownsamplingReadsIteratorTestDataProvider")
    public void runPerSampleDownsamplingReadsIteratorTest( PerSampleDownsamplingReadsIteratorTest test ) {
        logger.warn("Running test: " + test);

        Utils.resetRandomGenerator();
        test.run();
    }
}
