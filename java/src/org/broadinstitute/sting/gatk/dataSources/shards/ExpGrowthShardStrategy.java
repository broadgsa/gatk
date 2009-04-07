package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 8:23:19 PM
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
 * Class LinearShard
 * <p/>
 * A linear strategy, very very similar to adaptive
 */
public class ExpGrowthShardStrategy extends ShardStrategy {

    // fixed size
    private long baseSize = 100000;
    private long currentExp = 1;

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    ExpGrowthShardStrategy(SAMSequenceDictionary dic, long startSize) {
        super(dic);
        this.baseSize = startSize;
        currentExp = 1;
    }

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param strat the shatter to convert from
     */
    ExpGrowthShardStrategy(ShardStrategy strat) {
        super(strat);
        this.baseSize = strat.nextShardSize();
        currentExp = 1;
    }

    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
        baseSize = size;
        currentExp = 1;
    }

    /**
     * This is how the various shards strategies implements their approach
     *
     * @return the next shard size
     */
    protected long nextShardSize() {
        return (long) Math.floor(Math.pow((double) baseSize, (double) currentExp));
    }

}
