package org.broadinstitute.sting.gatk.iterators;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;



/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BoundedReadIteratorTest
 * <p/>
 * tests for the bounded read iterator.
 */
public class BoundedReadIteratorTest extends BaseTest {

    /** the file list and the fasta sequence */
    private List<File> fl;
    private FastaSequenceFile2 seq;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        fl = new ArrayList<File>();

        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq.getSequenceDictionary());
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testBounding() {
        logger.warn("Executing testBounding");
        // the sharding strat.
        ShardStrategy strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);
        int count = 0;


        // setup the test files
        fl.add(new File(seqLocation + "/dirseq/analysis/cancer_exome/twoflowcell_sams/TCGA-06-0188.aligned.duplicates_marked.bam"));
        Reads reads = new Reads(fl);

        // our target read
        final long boundedReadCount = 100;
        long shardReadCount = 0;

        try {
            SAMDataSource data = new SAMDataSource(reads,true);

            // make sure we have a shard
            if (!strat.hasNext()) {
                fail("Our shatter didn't give us a single piece, this is bad");
            }
            Shard sd = strat.next();


            StingSAMIterator datum = data.seek(sd);
            StingSAMIterator datum2 = data.seek(sd);

            // check the reads in the shard
            for (SAMRecord r : datum) {
                shardReadCount++;

            }

            // create the bounded iterator
            BoundedReadIterator iter = new BoundedReadIterator(StingSAMIteratorAdapter.adapt(reads,datum2), boundedReadCount);

            // now see how many reads are in the bounded iterator
            int readCount = 0;
            for (SAMRecord r : iter) {
                readCount++;

            }

            // close the iterators
            datum.close();
            datum2.close();

            // check to see that the sizes are the same
            assertEquals(boundedReadCount,readCount);
            assertTrue(readCount < shardReadCount);


        }
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a SimpleDataSourceLoadException");
        }


    }
}
