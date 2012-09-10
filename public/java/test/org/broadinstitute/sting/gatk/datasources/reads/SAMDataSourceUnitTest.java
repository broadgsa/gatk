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

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.downsampling.DownsampleType;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingMethod;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.ReadTransformer;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSingleSampleReadStream;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

/**
 * <p/>
 * Class SAMDataSourceUnitTest
 * <p/>
 * The test of the SAMBAM simple data source.
 */
public class SAMDataSourceUnitTest extends BaseTest {

    private List<SAMReaderID> readers;
    private IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;


    /***********************************
     * Tests for the fillShard() method
     ***********************************/

    /**
     * Tests to ensure that the fillShard() method does not place shard boundaries at inappropriate places,
     * such as within an alignment start position
     */
    private static class SAMDataSourceFillShardBoundaryTest extends TestDataProvider {
        private int numContigs;
        private int numStacksPerContig;
        private int stackSize;
        private int numUnmappedReads;
        private DownsamplingMethod downsamplingMethod;

        private SAMFileHeader header;

        public SAMDataSourceFillShardBoundaryTest( int numContigs,
                                                   int numStacksPerContig,
                                                   int stackSize,
                                                   int numUnmappedReads,
                                                   int downsamplingTargetCoverage ) {
            super(SAMDataSourceFillShardBoundaryTest.class);

            this.numContigs = numContigs;
            this.numStacksPerContig = numStacksPerContig;
            this.stackSize = stackSize;
            this.numUnmappedReads = numUnmappedReads;

            this.downsamplingMethod = new DownsamplingMethod(DownsampleType.BY_SAMPLE, downsamplingTargetCoverage, null, true);

            setName(String.format("%s: numContigs=%d numStacksPerContig=%d stackSize=%d numUnmappedReads=%d downsamplingTargetCoverage=%d",
                                  getClass().getSimpleName(), numContigs, numStacksPerContig, stackSize, numUnmappedReads, downsamplingTargetCoverage));
        }

        public void run() {
            SAMDataSource dataSource = new SAMDataSource(Arrays.asList(createTestBAM()),
                                                         new ThreadAllocation(),
                                                         null,
                                                         new GenomeLocParser(header.getSequenceDictionary()),
                                                         false,
                                                         SAMFileReader.ValidationStringency.SILENT,
                                                         null,
                                                         downsamplingMethod,
                                                         new ValidationExclusion(),
                                                         new ArrayList<ReadFilter>(),
                                                         false);

            Assert.assertTrue(dataSource.usingExpandedShards());

            Iterable<Shard> shardIterator = dataSource.createShardIteratorOverAllReads(new ReadShardBalancer());

            SAMRecord readAtEndOfLastShard = null;

            for ( Shard shard : shardIterator ) {
                int numContigsThisShard = 0;
                SAMRecord lastRead = null;

                for ( SAMRecord read : shard.iterator() ) {
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
        }

        private SAMReaderID createTestBAM() {
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

            File testBAMFile;
            try {
                testBAMFile = File.createTempFile("SAMDataSourceFillShardBoundaryTest", ".bam");
                testBAMFile.deleteOnExit();
            }
            catch ( IOException e ) {
                throw new ReviewedStingException(String.format("Failed to create temp bam file for test %s. %s", this, e.getMessage()));
            }

            SAMFileWriter bamWriter = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, testBAMFile);
            for ( SAMRecord read : artificialReads ) {
                bamWriter.addAlignment(read);
            }
            bamWriter.close();

            return new SAMReaderID(testBAMFile, new Tags());
        }
    }

    @DataProvider(name = "SAMDataSourceFillShardTestDataProvider")
    public Object[][] createSAMDataSourceFillShardBoundaryTests() {
        // Take downsampling out of the equation for these tests -- we are only interested in whether the
        // shard boundaries occur at the right places in the read stream, and removing downsampling as a
        // factor simplifies that task (note that we still need to provide a specific downsampling method with
        // experimental downsampling enabled to trigger the shard expansion behavior, for now)
        int downsamplingTargetCoverage = ReadShard.MAX_READS * 10;

        for ( int numContigs = 1; numContigs <= 3; numContigs++ ) {
            for ( int numStacksPerContig : Arrays.asList(1, 2, 4) ) {
                // Use crucial read shard boundary values as the stack sizes
                for ( int stackSize : Arrays.asList(ReadShard.MAX_READS / 2, ReadShard.MAX_READS / 2 + 10, ReadShard.MAX_READS, ReadShard.MAX_READS - 1, ReadShard.MAX_READS + 1, ReadShard.MAX_READS * 2) ) {
                    for ( int numUnmappedReads : Arrays.asList(0, ReadShard.MAX_READS / 2, ReadShard.MAX_READS * 2) ) {
                        new SAMDataSourceFillShardBoundaryTest(numContigs, numStacksPerContig, stackSize, numUnmappedReads, downsamplingTargetCoverage);
                    }
                }
            }
        }

        return SAMDataSourceFillShardBoundaryTest.getTests(SAMDataSourceFillShardBoundaryTest.class);
    }

    // TODO: re-enable these tests once the issues with filepointer ordering + the downsamplers are worked out
    @Test(dataProvider = "SAMDataSourceFillShardTestDataProvider", enabled = false)
    public void testSAMDataSourceFillShard( SAMDataSourceFillShardBoundaryTest test ) {
        logger.warn("Running test: " + test);

        test.run();
    }


    // TODO: the legacy tests below should really be replaced with a more comprehensive suite of tests for SAMDataSource

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeMethod
    public void doForEachTest() throws FileNotFoundException {
        readers = new ArrayList<SAMReaderID>();

        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(b36KGReference));
        genomeLocParser = new GenomeLocParser(seq.getSequenceDictionary());
    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @AfterMethod
    public void undoForEachTest() {
        seq = null;
        readers.clear();
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testLinearBreakIterateAll() {
        logger.warn("Executing testLinearBreakIterateAll");

        // setup the data
        readers.add(new SAMReaderID(new File(validationDataLocation+"/NA12878.chrom6.SLX.SRP000032.2009_06.selected.bam"),new Tags()));

        // the sharding strat.
        SAMDataSource data = new SAMDataSource(readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false);

        Iterable<Shard> strat = data.createShardIteratorOverMappedReads(seq.getSequenceDictionary(),new LocusShardBalancer());
        int count = 0;

        try {
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                GenomeLoc firstLocus = sh.getGenomeLocs().get(0), lastLocus = sh.getGenomeLocs().get(sh.getGenomeLocs().size()-1);
                logger.debug("Start : " + firstLocus.getStart() + " stop : " + lastLocus.getStop() + " contig " + firstLocus.getContig());
                logger.debug("count = " + count);
                StingSAMIterator datum = data.seek(sh);

                // for the first couple of shards make sure we can see the reads
                if (count < 5) {
                    for (SAMRecord r : datum) {
                    }
                    readCount++;
                }
                datum.close();

                // if we're over 100 shards, break out
                if (count > 100) {
                    break;
                }
            }
        }
        catch (UserException.CouldNotReadInputFile e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a UserException.CouldNotReadInputFile exception");
        }
    }

    /** Test that we clear program records when requested */
    @Test
    public void testRemoveProgramRecords() {
        logger.warn("Executing testRemoveProgramRecords");

        // setup the data
        readers.add(new SAMReaderID(new File(b37GoodBAM),new Tags()));

        // use defaults
        SAMDataSource data = new SAMDataSource(readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false);

        List<SAMProgramRecord> defaultProgramRecords = data.getHeader().getProgramRecords();
        assertTrue(defaultProgramRecords.size() != 0, "testRemoveProgramRecords: No program records found when using default constructor");

        boolean removeProgramRecords = false;
        data = new SAMDataSource(readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1,
                removeProgramRecords);

        List<SAMProgramRecord> dontRemoveProgramRecords = data.getHeader().getProgramRecords();
        assertEquals(dontRemoveProgramRecords, defaultProgramRecords, "testRemoveProgramRecords: default program records differ from removeProgramRecords = false");

        removeProgramRecords = true;
        data = new SAMDataSource(readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                SAMFileReader.ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1,
                removeProgramRecords);

        List<SAMProgramRecord> doRemoveProgramRecords = data.getHeader().getProgramRecords();
        assertTrue(doRemoveProgramRecords.isEmpty(), "testRemoveProgramRecords: program records not cleared when removeProgramRecords = true");
    }
}
