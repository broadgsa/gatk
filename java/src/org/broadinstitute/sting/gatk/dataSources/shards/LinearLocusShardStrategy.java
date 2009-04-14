package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

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
 * allows you to change the sharding length as you traverse
 */
class LinearLocusShardStrategy extends LocusShardStrategy {

    // default the next size to 100,000
    private long nextShardSize = 100000;

    /**
     * the constructor, taking a seq dictionary to parse out contigs
     *
     * @param dic the seq dictionary
     */
    LinearLocusShardStrategy(SAMSequenceDictionary dic, long startSize) {
        super(dic);
        this.nextShardSize = startSize;
    }

    /**
     * the constructor, constructing from another Locus strategy
     *
     * @param strat the shatter to convert from
     */
    LinearLocusShardStrategy(LocusShardStrategy strat) {
        super(strat);
        this.nextShardSize = strat.nextShardSize();
    }


    /**
     * The constructor, for a genomic list, start size, and a reference dictionary 
     * @param dic the reference dictionary
     * @param startSize the starting size of the shard
     * @param lst locations to iterate from
     */
    LinearLocusShardStrategy(SAMSequenceDictionary dic, long startSize, List<GenomeLoc> lst) {
        super(dic, lst);
        this.nextShardSize = startSize;
    }


    /**
     * set the next shards size
     *
     * @param size adjust the next size to this
     */
    public void adjustNextShardSize(long size) {
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
