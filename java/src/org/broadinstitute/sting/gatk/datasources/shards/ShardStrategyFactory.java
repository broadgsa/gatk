package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;

import java.io.File;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 7:09:22 PM
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
 * @date Apr 6, 2009
 * <p/>
 * Class ShardStrategyFactory
 * <p/>
 * The Shard Strategy Factory,  use this class to create and transfer shard strategies
 * between different approaches.
 */
public class ShardStrategyFactory {
    public enum SHATTER_STRATEGY {
        LINEAR,
        READS,
        INTERVAL,
        MONOLITHIC,   // Put all of the available data into one shard.
        LOCUS_EXPERIMENTAL,
        READS_EXPERIMENTAL
    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize) {
        return ShardStrategyFactory.shatter(readsDataSource, referenceDataSource, strat, dic, startingSize, -1L);
    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat               what's our strategy - SHATTER_STRATEGY type
     * @param dic                 the seq dictionary
     * @param startingSize        the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, long limitByCount) {
        switch (strat) {
            case LINEAR:
                return new LinearLocusShardStrategy(dic, startingSize, limitByCount);
            case READS:
                return new ReadDelimitedReadShardStrategy(startingSize, limitByCount);
            case INTERVAL:
                throw new StingException("Requested trategy: " + strat + " doesn't work with the limiting count (-M) command line option");
            case LOCUS_EXPERIMENTAL:
                return new IndexDelimitedLocusShardStrategy(readsDataSource,referenceDataSource,null);
            case READS_EXPERIMENTAL:
                return new BlockDelimitedReadShardStrategy(readsDataSource,null);
            default:
                throw new StingException("Strategy: " + strat + " isn't implemented for this type of shatter request");
        }

    }


    /**
     * get a new shatter strategy
     *
     * @param readsDataSource     File pointer to BAM.
     * @param referenceDataSource File pointer to reference.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return a shard strategy capable of dividing input data into shards.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocSortedSet lst) {
        return ShardStrategyFactory.shatter(readsDataSource, referenceDataSource, strat, dic, startingSize, lst, -1l);

    }

    /**
     * get a new shatter strategy
     *
     * @param readsDataSource The reads used to shatter this file.
     * @param referenceDataSource The reference used to shatter this file.
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return A strategy for shattering this data.
     */
    static public ShardStrategy shatter(SAMDataSource readsDataSource, IndexedFastaSequenceFile referenceDataSource, SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize, GenomeLocSortedSet lst, long limitDataCount) {
        switch (strat) {
            case LINEAR:
                return new LinearLocusShardStrategy(dic, startingSize, lst, limitDataCount);
            case INTERVAL:
                return new IntervalShardStrategy(startingSize, lst, Shard.ShardType.LOCUS_INTERVAL);
            case READS:
                return new IntervalShardStrategy(startingSize, lst, Shard.ShardType.READ_INTERVAL);
            case LOCUS_EXPERIMENTAL:
                return new IndexDelimitedLocusShardStrategy(readsDataSource,referenceDataSource,lst);
            case READS_EXPERIMENTAL:
                return new BlockDelimitedReadShardStrategy(readsDataSource,lst);
            default:
                throw new StingException("Strategy: " + strat + " isn't implemented");
        }

    }

}
