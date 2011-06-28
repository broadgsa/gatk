package org.broadinstitute.sting.gatk.datasources.reads;

import static org.testng.Assert.fail;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.datasources.reads.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.reads.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.reads.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
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
 * Class SAMBAMDataSourceUnitTest
 * <p/>
 * The test of the SAMBAM simple data source.
 */
public class SAMBAMDataSourceUnitTest extends BaseTest {

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
        seq = new CachingIndexedFastaSequenceFile(new File(hg18Reference));
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
        SAMDataSource data = new SAMDataSource(readers,genomeLocParser);
        ShardStrategy strat = ShardStrategyFactory.shatter(data,seq,ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL, seq.getSequenceDictionary(), 100000,genomeLocParser);
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


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testMergingTwoBAMFiles() {
        logger.warn("Executing testMergingTwoBAMFiles");

        // setup the test files
        readers.add(new SAMReaderID(new File(validationDataLocation + "/NA12878.chrom6.SLX.SRP000032.2009_06.selected.bam"),new Tags()));

        // the sharding strat.
        SAMDataSource data = new SAMDataSource(readers,genomeLocParser);
        ShardStrategy strat = ShardStrategyFactory.shatter(data,seq,ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL, seq.getSequenceDictionary(), 100000,genomeLocParser);

        ArrayList<Integer> readcountPerShard = new ArrayList<Integer>();
        ArrayList<Integer> readcountPerShard2 = new ArrayList<Integer>();

        // count up the first hundred shards
        int shardsToCount = 100;
        int count = 0;

        try {
            for (Shard sh : strat) {
                int readCount = 0;
                count++;
                if (count > shardsToCount) {
                    break;
                }

                StingSAMIterator datum = data.seek(sh);

                for (SAMRecord r : datum) {
                    readCount++;

                }
                readcountPerShard.add(readCount);
                logger.debug("read count = " + readCount);
                datum.close();
            }
        }
        catch (UserException.CouldNotReadInputFile e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a UserException.CouldNotReadInputFile exception");
        }


        // setup the data and the counter before our second run
        readers.clear();
        readers.add(new SAMReaderID(new File(validationDataLocation + "/NA12878.chrom6.SLX.SRP000032.2009_06.selected.bam"),new Tags()));
        readers.add(new SAMReaderID(new File(validationDataLocation + "/NA12878.chrom6.SLX.SRP000032.2009_06.selected.bam"),new Tags()));

        count = 0;
        // the sharding strat.
        data = new SAMDataSource(readers,genomeLocParser);
        strat = ShardStrategyFactory.shatter(data,seq,ShardStrategyFactory.SHATTER_STRATEGY.LOCUS_EXPERIMENTAL, seq.getSequenceDictionary(), 100000, genomeLocParser);

        logger.debug("Pile two:");
        try {
            for (Shard sh : strat) {
                int readCount = 0;
                count++;

                // can we leave?
                if (count > shardsToCount) {
                    break;
                }

                StingSAMIterator datum = data.seek(sh);

                for (SAMRecord r : datum) {
                    readCount++;
                }

                readcountPerShard2.add(readCount);
                logger.debug("read count = " + readCount);
                datum.close();
            }
        }
        catch (UserException.CouldNotReadInputFile e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should get a UserException.CouldNotReadInputFile exception");
        }

        /*int pos = 0;
        for (; pos < 100; pos++) {
            if (!readcountPerShard.get(pos).equals(readcountPerShard2.get(pos))) {
                fail("Shard number " + pos + " in the two approaches had different read counts, " + readcountPerShard.get(pos) + " and " + readcountPerShard2.get(pos));
            }
        } */

    }




}
