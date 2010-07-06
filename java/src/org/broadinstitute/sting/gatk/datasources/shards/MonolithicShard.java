package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Collections;
import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

/**
 * A single, monolithic shard bridging all available data.
 * @author mhanna
 * @version 0.1
 */
public class MonolithicShard implements Shard {
    /**
     * What type of MonolithicShard is this?  Read or locus?
     */
    private final ShardType shardType;

    /**
     * Locations.  For the monolithic shard, should be a list of all available contigs in the reference.
     */
    private final List<GenomeLoc> locs = new ArrayList<GenomeLoc>();

    /**
     * Creates a new monolithic shard of the given type.
     * @param shardType Type of the shard.  Must be either read or locus; cannot be intervalic.
     * @param sequenceDictionary the sequence dictionary from which to derive contig info. 
     */
    public MonolithicShard(ShardType shardType, SAMSequenceDictionary sequenceDictionary) {
        if(shardType != ShardType.LOCUS && shardType != ShardType.READ)
            throw new StingException("Invalid shard type for monolithic shard: " + shardType);
        this.shardType = shardType;
        for(SAMSequenceRecord sequenceRecord: sequenceDictionary.getSequences())
            locs.add(GenomeLocParser.createGenomeLoc(sequenceRecord.getSequenceName(),1,sequenceRecord.getSequenceLength()));
    }

    /**
     * Returns null, indicating that (in this case) the entire genome is covered.
     * @return null.
     */
    public List<GenomeLoc> getGenomeLocs() {
        return locs;
    }

    /**
     * Reports the type of monolithic shard.
     * @return Type of monolithic shard.
     */
    public ShardType getShardType() {
        return shardType;
    }

    /**
     * String representation of this shard.
     * @return "entire genome".
     */
    @Override
    public String toString() {
        return "entire genome";    
    }
}
