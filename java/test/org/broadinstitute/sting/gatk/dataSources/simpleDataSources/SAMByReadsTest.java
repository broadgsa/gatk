package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
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
    private List<String> fl;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        fl = new ArrayList<String>();

        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testTotalReadCount() {
        logger.warn("Executing testTotalReadCount");
        // the sharding strat.
        //ShardStrategy strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);

        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        // setup the test files
        fl.add("/humgen/gsa-scr1/GATK_Data/Validation_Data/10035.5.clean.bam");

        // make sure we don't see dupes
        HashSet<String> names = new HashSet<String>();

        
        try {
            SAMDataSource data = new SAMDataSource(fl);
            final int targetReadCount = 500000;

            // check the total read count
            final int totalReads = 1782980;

            int readsSeen = 0;
            BoundedReadIterator iter;


            while ((iter = data.seek(targetReadCount)) != null) {
                int readcnt = 0;
              
                for (SAMRecord r : iter) {
                    String readName = r.getReadName();
                    if (names.contains(readName)) {
                        fail("We saw read " + readName + " twice");
                    }
                    names.add(readName);
                    readcnt++;
                }

                
                readsSeen += readcnt;
                logger.warn("Seen " + readsSeen + " reads.");

            }
            // make sure we've seen all the reads
            assertEquals(totalReads,readsSeen);
        }
        
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }


    }
}
