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
 * Class ShardBuilder
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class ShardStrategyFactory {
    public enum SHATTER_STRATEGY {
        ADAPTIVE, LINEAR
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
        ShardStrategy d = null;
        switch (strat) {
            case ADAPTIVE:
                d = new AdaptiveShardStrategy(dic, startingSize);
            default:
                d = new LinearShardStrategy(dic, startingSize); // default
        }
        return d;
    }

    /**
     * if you know what you want
     *
     * @param dic          the seq dictionary
     * @param startingSize the starting size
     * @return
     */
    static public AdaptiveShardStrategy getAdaptiveShard(SAMSequenceDictionary dic, long startingSize) {
        return new AdaptiveShardStrategy(dic, startingSize);
    }

}
