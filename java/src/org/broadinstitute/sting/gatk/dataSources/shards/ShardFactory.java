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
public class ShardFactory {
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
    static public Shard shatter(SHATTER_STRATEGY strat, SAMSequenceDictionary dic, long startingSize) {
        Shard d = null;
        switch (strat) {
            case ADAPTIVE:
                d = new AdaptiveShard(dic, startingSize);
            default:
                d = new LinearShard(dic, startingSize); // default
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
    static public AdaptiveShard getAdaptiveShard(SAMSequenceDictionary dic, long startingSize) {
        return new AdaptiveShard(dic, startingSize);
    }

}
