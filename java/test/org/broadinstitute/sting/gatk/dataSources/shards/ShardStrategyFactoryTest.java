package org.broadinstitute.sting.gatk.dataSources.shards;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.FastaSequenceFile2;
import org.junit.*;

import java.io.File;
import java.util.ArrayList;

/**
 *
 * User: aaron
 * Date: Apr 8, 2009
 * Time: 11:31:04 AM
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
 * Class ShardFactoryTest
 * <p/>
 * Tests the shard strategy factory.  This tests the whole sharding interface, and should be
 * split in the future into seperate test cases.
 * TODO: split out for the seperate sharding classes
 */
public class ShardStrategyFactoryTest extends BaseTest {

    private static FastaSequenceFile2 seq;

    /**
     * This function (because of the @BeforeClass tag) gets called only once ever,
     * before any tests are run
     */
    @BeforeClass
    public static void doBeforeAnyTests() {
        seq = new FastaSequenceFile2(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @AfterClass
    public static void doAfterAllTests() {

    }

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {

    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @After
    public void undoForEachTest() {

    }

    /** Tests that we got a string parameter in correctly */
    @Test
    public void testFullGenomeCycle() {
        GenomeLoc.setupRefContigOrdering(seq.getSequenceDictionary());

        ShardStrategy strategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 100000);
        int shardCount = 0;
        try {

            for (Shard s : strategy) {
                GenomeLoc l = s.getGenomeLoc();
                //logger.debug("Shard start: " + l.getStart() + " stop " + l.getStop() + " contig " + l.getContig());
                shardCount++;
            }

            // check to make sure we got apple shards
            //logger.debug("shardCount : " + shardCount + " seq size = " + seq.getSequenceDictionary().size());

        } catch (Exception e) {
            e.printStackTrace();
            fail("We Shouldn't of seen an exception! : " + e.getMessage() + "; shard count " + shardCount);
        }
    }


    /** Tests that we got a string parameter in correctly */
    @Test
    public void testIntervalGenomeCycle() throws InterruptedException {
        SAMSequenceDictionary dic = seq.getSequenceDictionary();
        SAMSequenceRecord s = dic.getSequence(1);
        // Character stream writing

       
        int stop = s.getSequenceLength();
        int size = 10000;
        int location = 1;
        GenomeLoc.setupRefContigOrdering(dic);
        logger.debug("done to sleep");
        // keep track of the number of genome locs we build
        int genomeLocs = 0;
        ArrayList<GenomeLoc> locations = new ArrayList<GenomeLoc>();
        logger.debug("done to sleep2");
        try {
            while (location + size < stop) {
            logger.debug("s = " + s.getSequenceName() + " " + location + " " + size);
            // lets make up some fake locations
            GenomeLoc gl = new GenomeLoc(s.getSequenceName(), location, location + size - 1);
            logger.debug("loc = " + location);

            // let's move the location up, with a size space
            location += (size * 2);

            // add our current location to the list
            locations.add(gl);

            // add another genome location
            ++genomeLocs;
        }
        } catch (Exception e) {
            e.printStackTrace();
        }
        logger.debug("Location count = " + genomeLocs);
        ShardStrategy strategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 5000, locations);
        int shardCount = 0;
        try {
            for (Shard sh : strategy) {
                GenomeLoc l = sh.getGenomeLoc();

                logger.debug("Shard start: " + l.getStart() + " stop " + l.getStop() + " contig " + l.getContig());
                shardCount++;
            }

             logger.debug("Shard count = " + shardCount); 
            assertEquals(shardCount, genomeLocs * 2);

        } catch (Exception e) {
            e.printStackTrace();
            fail("testIntervalGenomeCycle: ne exception expected");
        }
    }

}
