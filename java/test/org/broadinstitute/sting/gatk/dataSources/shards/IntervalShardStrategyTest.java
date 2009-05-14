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
 * Date: May 14, 2009
 * Time: 3:52:57 PM
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
 * @date May 14, 2009
 * <p/>
 * Class LocusWindowShardStrategyTest
 * <p/>
 * LocusWindowShardStrategy tests
 */
public class IntervalShardStrategyTest extends BaseTest {

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
    public void testIntervalGenomeCycle() throws InterruptedException {
        logger.warn("Executing testIntervalGenomeCycle");

        SAMSequenceDictionary dic = seq.getSequenceDictionary();


        // setup a list of genome locs that represent the whole file
        SAMSequenceRecord s = dic.getSequence(1);
        int stop = s.getSequenceLength();
        int size = 10000;
        int location = 1;

        GenomeLoc.setupRefContigOrdering(dic);
        // keep track of the number of genome locs we build
        int genomeLocs = 0;
        ArrayList<GenomeLoc> locations = new ArrayList<GenomeLoc>();
        try {
            while (location + size < stop) {
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
        ShardStrategy strategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.INTERVAL, seq.getSequenceDictionary(), 0, locations);
        int shardCount = 0;
        try {
            for (Shard sh : strategy) {
                GenomeLoc l = sh.getGenomeLoc();
                GenomeLoc truth = locations.get(shardCount);
                if (l.compareTo(truth) != 0) {
                    String truthStr = truth.getContig() + ":" + truth.getStart() + ":" + truth.getStop();
                    String lStr = l.getContig() + ":" + l.getStart() + ":" + l.getStop();
                    fail("Genome loc " + truthStr + " doesn't equal " + lStr);
                }
                shardCount++;
            }
            assertEquals(shardCount, genomeLocs);

        } catch (Exception e) {
            e.printStackTrace();
            fail("testIntervalGenomeCycle: ne exception expected");
        }
    }

}
