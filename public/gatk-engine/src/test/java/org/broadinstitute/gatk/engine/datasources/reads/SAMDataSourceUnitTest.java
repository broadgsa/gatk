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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.commandline.Tags;
import org.broadinstitute.gatk.utils.ValidationExclusion;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.iterators.ReadTransformer;
import org.broadinstitute.gatk.utils.iterators.GATKSAMIterator;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
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

    // TODO: These legacy tests should really be replaced with a more comprehensive suite of tests for SAMDataSource

    private List<SAMReaderID> readers;
    private File referenceFile;
    private IndexedFastaSequenceFile seq;
    private GenomeLocParser genomeLocParser;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @BeforeMethod
    public void doForEachTest() throws FileNotFoundException {
        readers = new ArrayList<SAMReaderID>();

        // sequence
        referenceFile = new File(b36KGReference);
        seq = new CachingIndexedFastaSequenceFile(referenceFile);
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
        SAMDataSource data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false);

        Iterable<Shard> strat = data.createShardIteratorOverMappedReads(new LocusShardBalancer());
        int count = 0;

        try {
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                GenomeLoc firstLocus = sh.getGenomeLocs().get(0), lastLocus = sh.getGenomeLocs().get(sh.getGenomeLocs().size()-1);
                logger.debug("Start : " + firstLocus.getStart() + " stop : " + lastLocus.getStop() + " contig " + firstLocus.getContig());
                logger.debug("count = " + count);
                GATKSAMIterator datum = data.seek(sh);

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
        SAMDataSource data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false);

        List<SAMProgramRecord> defaultProgramRecords = data.getHeader().getProgramRecords();
        assertTrue(defaultProgramRecords.size() != 0, "testRemoveProgramRecords: No program records found when using default constructor");

        boolean removeProgramRecords = false;
        data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1,
                removeProgramRecords,
                false,
                null, IntervalMergingRule.ALL);

        List<SAMProgramRecord> dontRemoveProgramRecords = data.getHeader().getProgramRecords();
        assertEquals(dontRemoveProgramRecords, defaultProgramRecords, "testRemoveProgramRecords: default program records differ from removeProgramRecords = false");

        removeProgramRecords = true;
        data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1,
                removeProgramRecords,
                false,
                null, IntervalMergingRule.ALL);

        List<SAMProgramRecord> doRemoveProgramRecords = data.getHeader().getProgramRecords();
        assertTrue(doRemoveProgramRecords.isEmpty(), "testRemoveProgramRecords: program records not cleared when removeProgramRecords = true");
    }

    @Test(expectedExceptions = UserException.class)
    public void testFailOnReducedReads() {
        readers.add(new SAMReaderID(new File(privateTestDir + "old.reduced.bam"), new Tags()));

        SAMDataSource data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                false);
    }

    @Test(expectedExceptions = UserException.class)
    public void testFailOnReducedReadsRemovingProgramRecords() {
        readers.add(new SAMReaderID(new File(privateTestDir + "old.reduced.bam"), new Tags()));

        SAMDataSource data = new SAMDataSource(
                referenceFile,
                readers,
                new ThreadAllocation(),
                null,
                genomeLocParser,
                false,
                ValidationStringency.SILENT,
                null,
                null,
                new ValidationExclusion(),
                new ArrayList<ReadFilter>(),
                Collections.<ReadTransformer>emptyList(),
                false,
                (byte) -1,
                true,
                false,
                null, IntervalMergingRule.ALL);
    }
}
