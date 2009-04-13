package org.broadinstitute.sting.gatk.dataSources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 *
 * User: aaron
 * Date: Apr 7, 2009
 * Time: 1:19:49 PM
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
 * @date Apr 7, 2009
 * <p/>
 * Class Shard
 * <p/>
 * This is the base class for shards.  Right now it does little more then
 * wrap GenomeLoc (actually nothing more), but it's good to have the class
 * in place so it's easier to change guts later.
 */
public class LocusShard implements Shard {

    // currently our location
    final GenomeLoc mLoc;

    public LocusShard(GenomeLoc loc) {
        this.mLoc = loc;
    }

    /** @return the genome location represented by this shard */
    public GenomeLoc getGenomeLoc() {
        return mLoc;
    }

    /**
     * what kind of shard do we return
     *
     * @return ShardType, indicating the type
     */
    public ShardType getShardType() {
        return ShardType.LOCUS;
    }

    /**
     * return a shard representing the passed in GenomeLoc
     *
     * @return
     */
    public static LocusShard toShard(GenomeLoc loc) {
        return new LocusShard(loc);
    }
}
