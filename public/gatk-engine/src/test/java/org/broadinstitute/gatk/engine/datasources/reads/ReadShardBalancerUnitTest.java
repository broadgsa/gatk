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

import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.ArtificialSingleSampleReadStream;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class ReadShardBalancerUnitTest extends BaseTest {

    /**
     * Tests to ensure that ReadShardBalancer works as expected and does not place shard boundaries
     * at inappropriate places, such as within an alignment start position
     */
    private static class ReadShardBalancerTest extends TestDataProvider {
        private int numContigs;
        private int numStacksPerContig;
        private int stackSize;
        private int numUnmappedReads;
        private DownsamplingMethod downsamplingMethod;
        private int expectedReadCount;

        private SAMFileHeader header;
        private SAMReaderID testBAM;

        public ReadShardBalancerTest( int numContigs,
                                      int numStacksPerContig,
                                      int stackSize,
                                      int numUnmappedReads,
                                      int downsamplingTargetCoverage ) {
            super(ReadShardBalancerTest.class);

            this.numContigs = numContigs;
            this.numStacksPerContig = numStacksPerContig;
            this.stackSize = stackSize;
            this.numUnmappedReads = numUnmappedReads;

            this.downsamplingMethod = new DownsamplingMethod(DownsampleType.BY_SAMPLE, downsamplingTargetCoverage, null);
            this.expectedReadCount = Math.min(stackSize, downsamplingTargetCoverage) * numStacksPerContig * numContigs + numUnmappedReads;

            setName(String.format("%s: numContigs=%d numStacksPerContig=%d stackSize=%d numUnmappedReads=%d downsamplingTargetCoverage=%d",
                                  getClass().getSimpleName(), numContigs, numStacksPerContig, stackSize, numUnmappedReads, downsamplingTargetCoverage));
        }

        public void run() {
            createTestBAM();

            SAMDataSource dataSource = new SAMDataSource(null, // Reference not used in this test.
                                                         Arrays.asList(testBAM),
                                                         new ThreadAllocation(),
                                                         null,
                                                         new GenomeLocParser(header.getSequenceDictionary()),
                                                         false,
                                                         ValidationStringency.SILENT,
                                                         ReadShard.DEFAULT_MAX_READS,  // reset ReadShard.MAX_READS to ReadShard.DEFAULT_MAX_READS for each test
                                                         downsamplingMethod,
                                                         new ValidationExclusion(),
                                                         new ArrayList<ReadFilter>(),
                                                         false);

            Iterable<Shard> shardIterator = dataSource.createShardIteratorOverAllReads(new ReadShardBalancer());

            SAMRecord readAtEndOfLastShard = null;
            int totalReadsSeen = 0;

            for ( Shard shard : shardIterator ) {
                int numContigsThisShard = 0;
                SAMRecord lastRead = null;

                for ( SAMRecord read : shard.iterator() ) {
                    totalReadsSeen++;

                    if ( lastRead == null ) {
                        numContigsThisShard = 1;
                    }
                    else if ( ! read.getReadUnmappedFlag() && ! lastRead.getReferenceIndex().equals(read.getReferenceIndex()) ) {
                        numContigsThisShard++;
                    }

                    // If the last read from the previous shard is not unmapped, we have to make sure
                    // that no reads in this shard start at the same position
                    if ( readAtEndOfLastShard != null && ! readAtEndOfLastShard.getReadUnmappedFlag() ) {
                        Assert.assertFalse(readAtEndOfLastShard.getReferenceIndex().equals(read.getReferenceIndex()) &&
                                           readAtEndOfLastShard.getAlignmentStart() == read.getAlignmentStart(),
                                           String.format("Reads from alignment start position %d:%d are split across multiple shards",
                                                         read.getReferenceIndex(), read.getAlignmentStart()));
                    }

                    lastRead = read;
                }

                // There should never be reads from more than 1 contig in a shard (ignoring unmapped reads)
                Assert.assertTrue(numContigsThisShard == 1, "found a shard with reads from multiple contigs");

                readAtEndOfLastShard = lastRead;
            }

            Assert.assertEquals(totalReadsSeen, expectedReadCount, "did not encounter the expected number of reads");
        }

        private void createTestBAM() {
            header = ArtificialSAMUtils.createArtificialSamHeader(numContigs, 1, 100000);
            SAMReadGroupRecord readGroup = new SAMReadGroupRecord("foo");
            readGroup.setSample("testSample");
            header.addReadGroup(readGroup);
            ArtificialSingleSampleReadStream artificialReads = new ArtificialSingleSampleReadStream(header,
                                                                                                    "foo",
                                                                                                    numContigs,
                                                                                                    numStacksPerContig,
                                                                                                    stackSize,
                                                                                                    stackSize,
                                                                                                    1,
                                                                                                    100,
                                                                                                    50,
                                                                                                    150,
                                                                                                    numUnmappedReads);

            final File testBAMFile = createTempFile("SAMDataSourceFillShardBoundaryTest", ".bam");

            SAMFileWriter bamWriter = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, testBAMFile);
            for ( SAMRecord read : artificialReads ) {
                bamWriter.addAlignment(read);
            }
            bamWriter.close();

            testBAM =  new SAMReaderID(testBAMFile, new Tags());

            new File(testBAM.getSamFilePath().replace(".bam", ".bai")).deleteOnExit();
            new File(testBAM.getSamFilePath() + ".bai").deleteOnExit();
        }
    }

    @DataProvider(name = "ReadShardBalancerTestDataProvider")
    public Object[][] createReadShardBalancerTests() {
        for ( int numContigs = 1; numContigs <= 3; numContigs++ ) {
            for ( int numStacksPerContig : Arrays.asList(1, 2, 4) ) {
                // Use crucial read shard boundary values as the stack sizes
                for ( int stackSize : Arrays.asList(ReadShard.DEFAULT_MAX_READS / 2, ReadShard.DEFAULT_MAX_READS / 2 + 10, ReadShard.DEFAULT_MAX_READS, ReadShard.DEFAULT_MAX_READS - 1, ReadShard.DEFAULT_MAX_READS + 1, ReadShard.DEFAULT_MAX_READS * 2) ) {
                    for ( int numUnmappedReads : Arrays.asList(0, ReadShard.DEFAULT_MAX_READS / 2, ReadShard.DEFAULT_MAX_READS * 2) ) {
                        // The first value will result in no downsampling at all, the others in some downsampling
                        for ( int downsamplingTargetCoverage : Arrays.asList(ReadShard.DEFAULT_MAX_READS * 10, ReadShard.DEFAULT_MAX_READS, ReadShard.DEFAULT_MAX_READS / 2) ) {
                            new ReadShardBalancerTest(numContigs, numStacksPerContig, stackSize, numUnmappedReads, downsamplingTargetCoverage);
                        }
                    }
                }
            }
        }

        return ReadShardBalancerTest.getTests(ReadShardBalancerTest.class);
    }

    @Test(dataProvider = "ReadShardBalancerTestDataProvider")
    public void runReadShardBalancerTest( ReadShardBalancerTest test ) {
        logger.warn("Running test: " + test);

        test.run();
    }
}
