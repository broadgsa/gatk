package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.dataSources.shards.Shard;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.MergingSamRecordIterator2;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
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
 * The test of the SAMBAM simple data source.
 */
public class SAMBAMDataSourceTest extends BaseTest {

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

        // sequence
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());
    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @After
    public void undoForEachTest() {
        seq = null;
        fl.clear();
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testLinearBreakIterateAll() {
        logger.warn("Executing testLinearBreakIterateAll");
        // the sharding strat.
        ShardStrategy strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);
        int count = 0;

        // setup the data
        fl.add(oneKGLocation + "/pilot3/sams/NA12892.bam");


        try {
            SAMDataSource data = new SAMDataSource(fl);
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                logger.debug("Start : " + sh.getGenomeLoc().getStart() + " stop : " + sh.getGenomeLoc().getStop() + " contig " + sh.getGenomeLoc().getContig());
                logger.debug("count = " + count);
                MergingSamRecordIterator2 datum = data.seek(sh.getGenomeLoc());

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
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a SimpleDataSourceLoadException");
        }
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testMergingTwoBAMFiles() {
        logger.warn("Executing testMergingTwoBAMFiles");
        // the sharding strat.
        ShardStrategy strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);


        // setup the test files
        fl.add(seqLocation + "/dirseq/analysis/cancer_exome/twoflowcell_sams/TCGA-06-0188.aligned.duplicates_marked.bam");

        ArrayList<Integer> readcountPerShard = new ArrayList<Integer>();
        ArrayList<Integer> readcountPerShard2 = new ArrayList<Integer>();

        // count up the first hundred shards
        int shardsToCount = 100;
        int count = 0;

        try {
            SAMDataSource data = new SAMDataSource(fl);
            for (Shard sh : strat) {
                int readCount = 0;
                count++;
                if (count > shardsToCount) {
                    break;
                }

                MergingSamRecordIterator2 datum = data.seek(sh.getGenomeLoc());

                for (SAMRecord r : datum) {
                    readCount++;

                }
                readcountPerShard.add(readCount);
                logger.debug("read count = " + readCount);
                datum.close();
            }
        }
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a SimpleDataSourceLoadException");
        }


        // setup the data and the counter before our second run
        fl.clear();
        fl.add(seqLocation + "/dirseq/analysis/cancer_exome/twoflowcell_sams/TCGA-06-0188-01A-01W.aligned.duplicates_marked.bam");
        fl.add(seqLocation + "/dirseq/analysis/cancer_exome/twoflowcell_sams/TCGA-06-0188-10B-01W.aligned.duplicates_marked.bam");
        count = 0;
        // the sharding strat.
        strat = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);

        logger.debug("Pile two:");
        try {
            SAMDataSource data = new SAMDataSource(fl);
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                // can we leave?
                if (count > shardsToCount) {
                    break;
                }

                MergingSamRecordIterator2 datum = data.seek(sh.getGenomeLoc());

                for (SAMRecord r : datum) {
                    readCount++;
                }

                readcountPerShard2.add(readCount);
                logger.debug("read count = " + readCount);
                datum.close();
            }
        }
        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a SimpleDataSourceLoadException");
        }

        /*int pos = 0;
        for (; pos < 100; pos++) {
            if (!readcountPerShard.get(pos).equals(readcountPerShard2.get(pos))) {
                fail("Shard number " + pos + " in the two approaches had different read counts, " + readcountPerShard.get(pos) + " and " + readcountPerShard2.get(pos));
            }
        } */

    }




}
