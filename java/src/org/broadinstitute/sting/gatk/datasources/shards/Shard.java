package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.ReadMetrics;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.HasGenomeLocation;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.Serializable;
import java.util.List;
/**
 *
 * User: aaron
 * Date: Apr 10, 2009
 * Time: 5:00:27 PM
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
 * @date Apr 10, 2009
 * <p/>
 * Interface Shard
 * <p/>
 * The base abstract class for shards.
 */
public abstract class Shard implements HasGenomeLocation {
    public enum ShardType {
        READ, LOCUS
    }

    protected final GenomeLocParser parser; // incredibly annoying!

    /**
     * What type of MonolithicShard is this?  Read or locus?
     */
    protected final ShardType shardType;

    /**
     * Locations.  For the monolithic shard, should be a list of all available contigs in the reference.
     */
    protected final List<GenomeLoc> locs;

    /**
     * Statistics about which reads in this shards were used and which were filtered away.
     */
    protected final ReadMetrics readMetrics = new ReadMetrics();

    public Shard(GenomeLocParser parser, ShardType shardType, List<GenomeLoc> locs) {
        this.locs = locs;
        this.parser = parser;
        this.shardType = shardType;
    }

    /**
     * If isUnmapped is true, than getGenomeLocs by
     * definition will return a singleton list with a GenomeLoc.UNMAPPED
     *
     * Can return null, indicating that the entire genome is covered.
     *
     * @return the genome location represented by this shard
     */
    public List<GenomeLoc> getGenomeLocs() {
        return locs;
    }

    /**
     * Returns the span of the genomeLocs comprising this shard
     * @param
     * @return
     */
    public GenomeLoc getLocation() {
        if ( getGenomeLocs() == null )
            return GenomeLoc.WHOLE_GENOME;

        int start = Integer.MAX_VALUE;
        int stop = Integer.MIN_VALUE;
        String contig = null;

        for ( GenomeLoc loc : getGenomeLocs() ) {
            if ( GenomeLoc.isUnmapped(loc) )
                // special case the unmapped region marker, just abort out
                return loc;
            contig = loc.getContig();
            if ( loc.getStart() < start ) start = loc.getStart();
            if ( loc.getStop() > stop ) stop = loc.getStop();
        }

        return parser.createGenomeLoc(contig, start, stop);
    }


    /**
     * what kind of shard do we return
     * @return ShardType, indicating the type
     */
    public ShardType getShardType() {
        return shardType;
    }

    /**
     * Does any releasing / aggregation required when the shard is through being processed.
     */
    public void close() {
        ; // by default don't do anything
    }

    /**
     * Gets required configuration for validating and filtering reads.
     * @return read configuration properties.
     */
    public abstract ReadProperties getReadProperties();

    /**
     * Gets the runtime metrics associated with this shard.
     * Retrieves a storage space of metrics about number of reads included, filtered, etc.
     * @return Storage space for metrics.
     */
    public ReadMetrics getReadMetrics() {
        return readMetrics;
    }
}
