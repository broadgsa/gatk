package org.broadinstitute.sting.gatk.iterators;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.dataSources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 14, 2009
 * Time: 5:48:48 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BoundedReadIteratorTest
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
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
        GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());
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
            SAMDataSource data = new SAMDataSource(reads);

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
