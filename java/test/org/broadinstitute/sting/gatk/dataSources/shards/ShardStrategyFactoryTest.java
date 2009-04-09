package org.broadinstitute.sting.gatk.dataSources.shards;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.utils.FastaSequenceFile2;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.junit.*;

import java.io.File;
import java.io.FileWriter;
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
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class ShardStrategyFactoryTest {

    FastaSequenceFile2 seq = null;

    /**
     * This function (because of the @BeforeClass tag) gets called only once ever,
     * before any tests are run
     */
    @BeforeClass
    public static void doBeforeAnyTests() {

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

    /** Tests that we got a string parameter in correctly */
    @Test
    public void testFullGenomeCycle() {
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
    public void testIntervalGenomeCycle() {
        SAMSequenceDictionary dic = seq.getSequenceDictionary();
        SAMSequenceRecord s = dic.getSequence(1);
        // Character stream writing

        System.err.println("Trying to sleep");
        try {
            Thread.sleep(5000);
        } catch (InterruptedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        int stop = s.getSequenceLength();
        int size = 10000;
        int location = 1;
        System.err.println("done to sleep");
        // keep track of the number of genome locs we build
        int genomeLocs = 0;
        ArrayList<GenomeLoc> locations = new ArrayList<GenomeLoc>();
        while (location + size < stop) {
            // lets make up some fake locations
            GenomeLoc gl = new GenomeLoc(s.getSequenceName(), location, location + size - 1);

            // let's move the location up, with a size space
            location += (size * 2);

            // add our current location to the list
            locations.add(gl);

            // add another genome location
            ++genomeLocs;
        }

        ShardStrategy strategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR, seq.getSequenceDictionary(), 5000, locations);
        int shardCount = 0;
        try {
            FileWriter writer = new FileWriter("myfile.txt");
            for (Shard sh : strategy) {
                GenomeLoc l = sh.getGenomeLoc();

                writer.write("Shard start: " + l.getStart() + " stop " + l.getStop() + " contig " + l.getContig());
                //logger.debug("Shard start: " + l.getStart() + " stop " + l.getStop() + " contig " + l.getContig());
                shardCount++;
            }
            writer.close();
            assertEquals(shardCount, genomeLocs * 2);

        } catch (Exception e) {
            e.printStackTrace();
            fail("testIntervalGenomeCycle: ne exception expected");
        }
    }

}
