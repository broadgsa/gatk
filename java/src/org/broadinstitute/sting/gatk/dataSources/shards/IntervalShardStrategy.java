package org.broadinstitute.sting.gatk.dataSources.shards;

import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;

import java.util.Iterator;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
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
public class IntervalShardStrategy implements ShardStrategy {

    /** our storage of the genomic locations they'd like to shard over */
    private final GenomeLocSortedSet regions;

    /** their prefered size of the shard, we can modify this based on what we see in the shards */
    private long size;

    /**
     * change the recommended shard size for the next shard we generate.  The code will do it's
     * best to respect this value, but there are no guarantees.
     *
     * @param size the next recommended shard size.
     */
    public void adjustNextShardSize( long size ) {
        this.size = size;
    }

    /**
     * construct the shard strategy from a seq dictionary, a shard size, and and genomeLocs
     *
     * @param size
     * @param locations
     */
    IntervalShardStrategy( long size, GenomeLocSortedSet locations ) {
        if (locations == null || locations.isEmpty()) {
            throw new StingException("IntervalShardStrategy: genomic regions list is empty.");
        }
        this.regions = locations.clone();
        this.size = size;
    }

    /**
     * the default constructor
     *
     * @param dict the sequence dictionary to use
     */
    IntervalShardStrategy( SAMSequenceDictionary dict, GenomeLocSortedSet locations ) {
        if (locations == null || locations.isEmpty()) {
            throw new StingException("IntervalShardStrategy: genomic regions list is empty.");
        }
        this.regions = locations.clone();
        this.size = size;
    }

    /**
     * returns true if there are additional shards
     *
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return ( !regions.isEmpty() );
    }

    /**
     * gets the next Shard
     *
     * @return the next shard
     */
    public Shard next() {
        if (( this.regions == null ) || ( regions.isEmpty() )) {
            throw new StingException("IntervalShardStrategy: genomic regions list is empty in next() function.");
        }

        // get the first region in the list
        GenomeLoc loc = regions.iterator().next();

        regions.removeRegion(loc);
        return new IntervalShard(loc);

    }

    /** we don't support the remove command */
    public void remove() {
        throw new UnsupportedOperationException("ShardStrategies don't support remove()");
    }

    /**
     * makes the ReadIntervalShard iterable, i.e. usable in a for loop.
     *
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }
}
