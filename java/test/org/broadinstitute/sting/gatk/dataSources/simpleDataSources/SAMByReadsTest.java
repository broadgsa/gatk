package org.broadinstitute.sting.gatk.dataSources.simpleDataSources;

import static junit.framework.Assert.fail;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.dataSources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.iterators.BoundedReadIterator;
import org.broadinstitute.sting.gatk.iterators.*;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.ArtificialSAMQueryIterator;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Collections;

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
        fl.add(new File("/humgen/gsa-scr1/GATK_Data/Validation_Data/index_test.bam"));
        reads = new Reads(fl);


    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testToUnmappedReads() {
        ArtificialResourcePool gen = new ArtificialResourcePool(1,10,100,1000);
        GenomeLoc.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        try {
            int unmappedReadsSeen = 0;
            int iterations = 0;
            SAMDataSource data = new SAMDataSource(reads,true);
            for (int x = 0; x < 10; x++) {
                ++iterations;
                PeekingStingIterator iter = ArtificialSAMUtils.unmappedReadIterator(1, 100, 10, 1000);
                BoundedReadIterator ret = data.toUnmappedReads(100, iter);
                // count the reads we've gotten back
                if (ret == null) {
                    fail("On iteration " + iterations + " we were returned a null pointer, after seeing " + unmappedReadsSeen + " reads out of a 1000");
                }
                while (ret.hasNext()) {
                    ret.next();
                    unmappedReadsSeen++;
                }
            }
            assertEquals(1000,unmappedReadsSeen);
        }

        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }


    }

    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testShardingOfReadsSize14() {
        ArtificialResourcePool gen = new ArtificialResourcePool(1,10,100,1000);
        GenomeLoc.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        targetReadCount = 14;
        try {
            int iterations = 0;
            int readCount = 0;
            SAMDataSource data = new SAMDataSource(reads,true);


            data.setResourcePool(gen);
            shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS, gen.getHeader().getSequenceDictionary(), targetReadCount);
            while (shardStrategy.hasNext()) {


                BoundedReadIterator ret = (BoundedReadIterator)data.seek(shardStrategy.next());
                assertTrue(ret != null);
                while (ret.hasNext()) {
                    ret.next();
                    readCount++;
                }
                ret.close();
                iterations++;
            }

            // assert that we saw 2000 reads
            assertEquals(2000,readCount);

            /**
             * this next assertion is based on the following logic:
             * 14 reads per shard = 8 shards per each 100 read chromosome
             * 10 chromosomes = 8 * 10 = 80
             * 1000 unmapped reads / 14 = 72
             * 1 iteration at the end to know we're done
             * 80 + 72 + 1 = 153
             */
            assertEquals(153,iterations);

        }

        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }


    }

    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testShardingOfReadsSize25() {
        ArtificialResourcePool gen = new ArtificialResourcePool(1,10,100,1000);
        GenomeLoc.setupRefContigOrdering(gen.getHeader().getSequenceDictionary());
        targetReadCount = 25;
        try {
            int iterations = 0;
            int readCount = 0;
            SAMDataSource data = new SAMDataSource(reads,true);


            data.setResourcePool(gen);
            shardStrategy = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS, gen.getHeader().getSequenceDictionary(), targetReadCount);
            while (shardStrategy.hasNext()) {


                BoundedReadIterator ret = (BoundedReadIterator)data.seek(shardStrategy.next());
                assertTrue(ret != null);
                while (ret.hasNext()) {
                    ret.next();
                    readCount++;
                }
                ret.close();
                iterations++;
            }

            // assert that we saw 2000 reads
            assertEquals(2000,readCount);

            /**
             * this next assertion is based on the following logic:
             * 25 reads per shard = 5 shards (1 on the end to realize we're done)
             * 10 chromosomes = 5 * 10 = 50
             * 1000 unmapped reads / 25 = 40
             * 1 iteration at the end to know we're done
             * 50 + 40 + 1 = 91
             */
            assertEquals(91,iterations);

        }

        catch (SimpleDataSourceLoadException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail("testLinearBreakIterateAll: We Should not get a SimpleDataSourceLoadException");
        }


    }


}

/**
 * use this to inject into SAMDataSource for testing
 */
class ArtificialResourcePool extends SAMIteratorPool {
    // How strict should we be with SAM/BAM parsing?
    protected SAMFileReader.ValidationStringency strictness = SAMFileReader.ValidationStringency.SILENT;

    // the header
    private SAMFileHeader header;
    private final SAMFileHeader.SortOrder sortOrder = SAMFileHeader.SortOrder.coordinate;

    public ArtificialResourcePool( int startingChr, int endingChr, int readCount, int readSize) {
        super( new Reads(Collections.<File>emptyList()),true );
        header = ArtificialSAMUtils.createArtificialSamHeader(( endingChr - startingChr ) + 1, startingChr, readCount + readSize);

    }

    @Override
    public StingSAMIterator iterator( GenomeLoc loc ) {
        ArtificialSAMQueryIterator iter = ArtificialSAMUtils.queryReadIterator(1, 10, 100, 1000);
        if (loc != null) {
            iter.queryContained(loc.getContig(), (int)loc.getStart(), (int)loc.getStop());
        }
        return iter;
    }

    /**
     * get the merged header
     *
     * @return the merged header
     */
    public SAMFileHeader getHeader() {
       return this.header;
    }
}