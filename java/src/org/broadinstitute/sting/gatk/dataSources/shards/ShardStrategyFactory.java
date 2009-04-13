package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

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
        LINEAR, EXPONENTIAL, READS
    }

    /**
     * get a new shatter strategy
     *
     * @param strat        what's our strategy - SHATTER_STRATEGY type
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return
     */
    static public ShardStrategy shatter(SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize) {
        switch (strat) {
            case LINEAR:
                return new LinearLocusShardStrategy(dic, startingSize);
            case EXPONENTIAL:
                return new ExpGrowthLocusShardStrategy(dic, startingSize);
            default:
                throw new RuntimeException("Strategy: " + strat + " isn't implemented");
        }

    }

    /**
     * convert between types
     *
     * @param strat       the strategy
     * @param convertFrom convert from this strategy
     * @return
     */
    static public ShardStrategy transitionToShardStrategy(SHATTER_STRATEGY strat, LocusShardStrategy convertFrom) {
        switch (strat) {
            case LINEAR:
                return new LinearLocusShardStrategy(convertFrom);
            case EXPONENTIAL:
                return new ExpGrowthLocusShardStrategy(convertFrom);
            default:
                throw new RuntimeException("Strategy: " + strat + " isn't implemented");

        }
    }

    /**
     * convert between types
     *
     * @param readCount the number of reads to include in each shard
     * @return
     */
    static public ShardStrategy shatterByReadCount(long readCount) {
        return null;    
    }

}
