package org.broadinstitute.sting.gatk.datasources.shards;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.*;
import static org.junit.Assert.assertTrue;

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
 */
public class ShardStrategyFactoryTest extends BaseTest {

    private SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(NUMBER_OF_CHROMOSOMES, STARTING_CHROMOSOME, CHROMOSOME_SIZE);
    private static final int NUMBER_OF_CHROMOSOMES = 5;
    private static final int STARTING_CHROMOSOME = 1;
    private static final int CHROMOSOME_SIZE = 1000;
    private GenomeLocSortedSet set = null;


    @Before
    public void setup() {
        GenomeLocParser.setupRefContigOrdering(header.getSequenceDictionary());
        set = new GenomeLocSortedSet();
    }

    @Test
    public void testReadNonInterval() {
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS,header.getSequenceDictionary(),100);
        assertTrue(st instanceof ReadShardStrategy);
    }

    @Test
    public void testReadInterval() {
        GenomeLoc l = GenomeLocParser.createGenomeLoc(0,1,100);
        set.add(l);
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.READS,header.getSequenceDictionary(),100,set);
        assertTrue(st instanceof IntervalShardStrategy);
    }

    @Test
    public void testLinearNonInterval() {
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,header.getSequenceDictionary(),100);
        assertTrue(st instanceof LinearLocusShardStrategy);
    }

     @Test
    public void testExpNonInterval() {
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.EXPONENTIAL,header.getSequenceDictionary(),100);
        assertTrue(st instanceof ExpGrowthLocusShardStrategy);
    }

    @Test
    public void testExpInterval() {
        GenomeLoc l = GenomeLocParser.createGenomeLoc(0,1,100);
        set.add(l);
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.EXPONENTIAL,header.getSequenceDictionary(),100,set);
        assertTrue(st instanceof ExpGrowthLocusShardStrategy);
    }

    @Test
    public void testLinearInterval() {
        GenomeLoc l = GenomeLocParser.createGenomeLoc(0,1,100);
        set.add(l);
        ShardStrategy st = ShardStrategyFactory.shatter(ShardStrategyFactory.SHATTER_STRATEGY.LINEAR,header.getSequenceDictionary(),100,set);
        assertTrue(st instanceof LinearLocusShardStrategy);
    }
    
}
