package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.utils.FastaSequenceFile2;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * User: aaron
 * Date: Apr 8, 2009
 * Time: 8:14:23 PM
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
 * @date Apr 8, 2009
 * <p/>
 * Class SAMBAMDataSourceTest
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class SAMBAMDataSourceTest {

    private List<String> fl;
    private FastaSequenceFile2 seq;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        fl = new ArrayList<String>();
        fl.add("/humgen/gsa-scr1/aaron/stink/NA12892.bam");

        // sequence
        seq = new FastaSequenceFile2(new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @After
    public void undoForEachTest() {

    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testLinearBreakIterateAll() {
        // the sharding strat.
        ShardStrategy strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);
        int count = 0;

        try {
            SAMBAMDataSource data = new SAMBAMDataSource(fl);
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                System.out.println("Start : " + sh.getGenomeLoc().getStart() + " stop : " + sh.getGenomeLoc().getStop() + " contig " + sh.getGenomeLoc().getContig());
                System.out.println("count = " + count);
                MergingSamRecordIterator2 datum = data.seek(sh.getGenomeLoc());

                // for the first couple of shards make sure we can see the reads
                if (count < 5) {
                    for (SAMRecord r : datum) {
                    }
                    readCount++;
                }
                datum.close();
            }
        }
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a SimpleDataSourceLoadException");
        }
    }

}
