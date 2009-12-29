package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
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
 */
public class SAMByReadsTest extends BaseTest {


    private List<File> fl;
    ShardStrategy shardStrategy;
    Reads reads;
    private int targetReadCount = 14;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        fl = new ArrayList<File>();

        // sequence
        //seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly17/v0/Homo_sapiens_assembly17.fasta"));
        //GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());

        // setup the test files
        fl.add(new File(validationDataLocation + "index_test.bam"));
        reads = new Reads(fl);
    }


    /**
     * Test out that we can shard the file and iterate over every read
     */
    @Test
    public void testToUnmappedReads() {
        ArtificialResourcePool gen = new ArtificialResourcePool(createArtificialSamHeader(1, 1, 1, 50),
                ArtificialSAMUtils.mappedAndUnmappedReadIterator(1, 1, 1, 10));

        GenomeLocParser.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        try {
            int unmappedReadsSeen = 0;
            int iterations = 0;

            SAMDataSource data = new SAMDataSource(reads);
            data.setResourcePool(gen);
            ++iterations;
            StingSAMIterator ret = data.toUnmappedReads(100);
            // count the reads we've gotten back
            if (ret == null) {
                fail("On iteration " + iterations + " we were returned a null pointer, after seeing " + unmappedReadsSeen + " reads out of a 1000");
            }
            while (ret.hasNext()) {
                ret.next();
                unmappedReadsSeen++;
            }
            assertEquals(10, unmappedReadsSeen);
        } catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }
    }

    /**
     * Test out that we can shard the file and iterate over every read
     */
    @Test
    public void testShardingOfReadsEvenSplit() {
        ArtificialResourcePool gen = new ArtificialResourcePool(createArtificialSamHeader(1, 1, 10, 50),
                ArtificialSAMUtils.mappedAndUnmappedReadIterator(1, 1, 10, 10));
        GenomeLocParser.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        targetReadCount = 5;
        try {
            int readCount = 0;
            SAMDataSource data = new SAMDataSource(reads);

            data.setResourcePool(gen);
            shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS, gen.getHeader().getSequenceDictionary(), targetReadCount);
            while (shardStrategy.hasNext()) {
                StingSAMIterator ret = data.seek(shardStrategy.next());
                assertTrue(ret != null);
                while (ret.hasNext()) {
                    ret.next();
                    readCount++;
                }
                ret.close();
            }
            assertEquals(20, readCount);
        } catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }
    }

    /**
     * Test out that we can shard the file and iterate over every read
     */
    @Test
    public void testShardingOfReadsOddRemainder() {
        ArtificialResourcePool gen = new ArtificialResourcePool(createArtificialSamHeader(1, 1, 10, 100),
                ArtificialSAMUtils.mappedAndUnmappedReadIterator(1, 1, 10, 10));
        GenomeLocParser.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        targetReadCount = 3;
        try {
            int readCount = 0;
            SAMDataSource data = new SAMDataSource(reads);


            data.setResourcePool(gen);
            shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS, gen.getHeader().getSequenceDictionary(), targetReadCount);
            while (shardStrategy.hasNext()) {


                StingSAMIterator ret = data.seek(shardStrategy.next());
                assertTrue(ret != null);
                while (ret.hasNext()) {
                    ret.next();
                    readCount++;
                }
                ret.close();
            }

            // assert that we saw 2000 reads
            assertEquals(20, readCount);

        } catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }
    }

    private SAMFileHeader createArtificialSamHeader(int startingChr, int endingChr, int readCount, int readSize) {
        return ArtificialSAMUtils.createArtificialSamHeader((endingChr - startingChr) + 1,
                startingChr,
                readCount + readSize);
    }
}

