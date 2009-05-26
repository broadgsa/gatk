package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;

import java.util.Iterator;
import java.util.List;

/**
 *
 * User: aaron
 * Date: May 21, 2009
 * Time: 4:13:53 PM
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
 *         <p/>
 *         Class ReadByIntervalShardStrategy
 *         <p/>
 *         Impliments the sharding strategy for reads, given a list
 *         of genomic locations.  Shards returned will be bounded by the interval,
 *         but each provided interval may be split into a number of smaller regions.
 */
public class ReadIntervalShardStrategy implements ShardStrategy {

    /** our storage of the genomic locations they'd like to shard over */
    private final GenomeLocSortedSet regions;

    /** their prefered size of the shard, we can modify this based on what we see in the shards */
    private long size;

    /** the sequence dictionary we'll use to lookup the contigs */
    private final SAMSequenceDictionary dict;

    /**
     * change the recommended shard size for the next shard we generate.  The code will do it's
     * best to respect this value, but there are no guarantees.
     *
     * @param size the next recommended shard size.
     */
    public void adjustNextShardSize(long size) {
        this.size = size;
    }

    /**
     * the default constructor
     *
     * @param dict  the sequence dictionary to use
     * @param size the read count to iterate over
     */
    ReadIntervalShardStrategy(SAMSequenceDictionary dict, long size, GenomeLocSortedSet locations) {
        if (locations == null || locations.isEmpty()) {
            throw new StingException("ReadIntervalShardStrategy: genomic regions list is empty.");
        }
        this.regions = locations.clone();
        this.size = size;
        this.dict = dict;
    }

    /**
     * returns true if there are additional shards
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return (!regions.isEmpty());
    }

    /**
     * gets the next Shard
     * @return the next shard
     */
    public Shard next() {
        if ((this.regions == null) || (regions.isEmpty())) {
            throw new StingException("ReadIntervalShardStrategy: genomic regions list is empty in next() function.");
        }

        // get the first region in the list
        GenomeLoc loc = regions.iterator().next();

        if (loc.getStop() - loc.getStart() <= this.size) {
            regions.removeRegion(loc);
            return new IntervalReadShard(loc);
        } else {
            GenomeLoc subLoc = new GenomeLoc(loc.getContigIndex(),loc.getStart(),loc.getStart()+size-1);
            regions.removeRegion(subLoc);
            return new IntervalReadShard(subLoc);
        }
        
    }

    /**
     * we don't support the remove command
     */
    public void remove() {
        throw new UnsupportedOperationException("ShardStrategies don't support remove()");
    }

    /**
     * makes the ReadIntervalShard iterable, i.e. usable in a for loop.
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }
}
