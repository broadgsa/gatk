package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;

/**
 *
 * User: aaron
 * Date: Apr 6, 2009
 * Time: 7:18:19 PM
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
 * Class AdaptiveShard
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
class AdaptiveShard extends Shard {

    // default the next size to 100,000
    private long nextShardSize = 100000;

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    AdaptiveShard(SAMSequenceDictionary dic, long startSize) {
        super(dic);
        this.nextShardSize = startSize;
    }

    public void setNextShardSize(long size) {
        nextShardSize = size;
    }

    /**
     * This is how the various shards strategies implements their approach
     *
     * @return the next shard size
     */
    protected long nextShardSize() {
        return nextShardSize;
    }
}
