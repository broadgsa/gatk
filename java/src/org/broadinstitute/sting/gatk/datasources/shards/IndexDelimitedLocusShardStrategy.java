package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;

import java.util.*;

import net.sf.samtools.Chunk;

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
 * A sharding strategy for loci based on reading of the index.
 */
public class IndexDelimitedLocusShardStrategy implements ShardStrategy {

    /** our storage of the genomic locations they'd like to shard over */
    private final SortedMap<GenomeLoc,List<Chunk>> locations = new TreeMap<GenomeLoc,List<Chunk>>();

    /**
     * construct the shard strategy from a seq dictionary, a shard size, and and genomeLocs
     * @param dataSource Data source from which to load index data.
     * @param locations List of locations for which to load data.
     */
    IndexDelimitedLocusShardStrategy(SAMDataSource dataSource, GenomeLocSortedSet locations) {
        for(GenomeLoc location: locations)
            this.locations.put(location,((BlockDrivenSAMDataSource)dataSource).getOverlappingFilePointers(location));
    }

    /**
     * returns true if there are additional shards
     *
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return ( !locations.isEmpty() );
    }

    /**
     * gets the next Shard
     *
     * @return the next shard
     */
    public IndexDelimitedLocusShard next() {
        if (( this.locations == null ) || ( locations.isEmpty() )) {
            throw new StingException("IntervalShardStrategy: genomic regions list is empty in next() function.");
        }

        // get the first region in the list
        GenomeLoc loc = locations.firstKey();
        List<Chunk> filePointers = locations.get(loc);
        locations.remove(loc);
        
        return new IndexDelimitedLocusShard(GenomeLocSortedSet.createSetFromList(Arrays.asList(loc)),filePointers); 
    }

    /** we don't support the remove command */
    public void remove() {
        throw new UnsupportedOperationException("ShardStrategies don't support remove()");
    }

    /**
     * makes the IntervalShard iterable, i.e. usable in a for loop.
     *
     * @return
     */
    public Iterator<Shard> iterator() {
        return this;
    }
}