package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 15, 2009
 * Time: 7:12:59 PM
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
 * @date Apr 15, 2009
 * <p/>
 * Class SAMByReadsTest
 * <p/>
 * Test sam traversal by reads FIX ME
 */
public class SAMByReadsTest extends BaseTest {

    private FastaSequenceFile2 seq;
    private List<File> fl;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        fl = new ArrayList<File>();

        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly17/v0/Homo_sapiens_assembly17.fasta"));
        GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testTotalReadCount() {
        logger.warn("Executing testTotalReadCount");
       
        // setup the test files
        fl.add(new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/index_test.bam"));
        Reads reads = new Reads(fl);

        final int targetReadCount = 5000;
        
        ShardStrategy shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS,seq.getSequenceDictionary(),targetReadCount);
        
        try {
            SAMDataSource data = new SAMDataSource(reads);

            // check the total read count
            final int totalReads = 10000;

            int readsSeen = 0;
            BoundedReadIterator iter;


            for (Shard sd : shardStrategy) {
                int readcnt = 0;

                iter = (BoundedReadIterator)data.seek(sd);

                for (SAMRecord r : iter) {

                    readcnt++;

                }

                
                readsSeen += readcnt;
                //logger.warn("Seen " + readsSeen + " reads.");

            }
            // make sure we've seen all the reads
            assertEquals(totalReads,readsSeen);
            logger.warn("Success " + readsSeen + " equals target count of " + totalReads);

        }
        
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }


    }
}
