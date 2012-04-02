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

import static org.testng.Assert.fail;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;

import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author aaron
 * @version 1.0
 * @date Apr 8, 2009
 * <p/>
 * Class SAMDataSourceUnitTest
 * <p/>
 * The test of the SAMBAM simple data source.
 */
public class SAMDataSourceUnitTest extends BaseTest {

    private List<SAMReaderID> readers;
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
}
