package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.util.List;

/**
 * This is the base class for locus shards.  Right now it does little more then
 * wrap GenomeLoc (actually nothing more), but it's good to have the class
 * in place so it's easier to change guts later.
 * @author aaron
 * @version 1.0
 * @date Apr 7, 2009
 */
public class LocusShard implements Shard {
    // currently our location
    protected final List<GenomeLoc> loci;

    public LocusShard(List<GenomeLoc> loci) {
        this.loci = loci;
    }

    /** @return the genome location represented by this shard */
    public List<GenomeLoc> getGenomeLocs() {
        return loci;
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
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return Utils.join(";",loci);
    }
}
