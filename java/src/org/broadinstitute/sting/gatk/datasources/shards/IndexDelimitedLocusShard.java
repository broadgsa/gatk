package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import net.sf.samtools.Chunk;

import java.util.List;


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
 * A shard that's delimited based on the index rather than
 */
public class IndexDelimitedLocusShard implements Shard {

    /**
     * a collection of genomic locations to interate over
     */
    private final GenomeLocSortedSet intervals;

    /**
     * A list of the chunks associated with this shard.
     */
    private final List<Chunk> chunks;

    IndexDelimitedLocusShard(GenomeLocSortedSet intervals, List<Chunk> chunks) {
        this.intervals = intervals;
        this.chunks = chunks;
    }

    /**
     * The locations represented by this shard.
     * @return the genome location represented by this shard
     */
    public List<GenomeLoc> getGenomeLocs() {
        return intervals.toList();
    }

    /**
     * Gets the chunks associated with this locus shard.
     * @return A list of the chunks to use when retrieving locus data.
     */
    public List<Chunk> getChunks() {
        return chunks;
    }

    /**
     * returns the type of shard, LOCUS_INTERVAL.
     * @return LOCUS_INTERVAL, indicating the shard type
     */
    public ShardType getShardType() {
        return ShardType.LOCUS_INTERVAL;
    }

    /**
     * String representation of this shard.
     * @return A string representation of the boundaries of this shard.
     */
    @Override
    public String toString() {
        return intervals.toString();
    }
}