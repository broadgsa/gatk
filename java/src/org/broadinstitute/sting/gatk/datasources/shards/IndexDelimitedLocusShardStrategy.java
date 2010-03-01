package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;

import java.util.*;

import net.sf.samtools.Chunk;
import net.sf.samtools.SAMFileReader2;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

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
    /**
     * The data source to use when performing this sharding.
     */
    private final BlockDrivenSAMDataSource dataSource;

    /** our storage of the genomic locations they'd like to shard over */
    private final List<FilePointer> filePointers = new ArrayList<FilePointer>();

    /**
     * An iterator through the available file pointers.
     */
    private final Iterator<FilePointer> filePointerIterator;

    /**
     * construct the shard strategy from a seq dictionary, a shard size, and and genomeLocs
     * @param dataSource Data source from which to load index data.
     * @param locations List of locations for which to load data.
     */
    IndexDelimitedLocusShardStrategy(SAMDataSource dataSource, GenomeLocSortedSet locations) {
        if(!(dataSource instanceof BlockDrivenSAMDataSource))
            throw new StingException("Cannot power an IndexDelimitedLocusShardStrategy with this data source.");

        List<GenomeLoc> intervals;
        if(locations == null) {
            // If no locations were passed in, shard the entire BAM file.
            SAMFileHeader header = dataSource.getHeader();
            intervals = new ArrayList<GenomeLoc>();

            for(SAMSequenceRecord sequenceRecord: header.getSequenceDictionary().getSequences())
                intervals.add(GenomeLocParser.createGenomeLoc(sequenceRecord.getSequenceName(),1,sequenceRecord.getSequenceLength()));    
        }
        else
            intervals = locations.toList();


        this.dataSource = (BlockDrivenSAMDataSource)dataSource;
        filePointers.addAll(IntervalSharder.shardIntervals(this.dataSource,intervals,this.dataSource.getNumIndexLevels()-1));
        filePointerIterator = filePointers.iterator();
    }

    /**
     * returns true if there are additional shards
     *
     * @return false if we're done processing shards
     */
    public boolean hasNext() {
        return filePointerIterator.hasNext();
    }

    /**
     * gets the next Shard
     *
     * @return the next shard
     */
    public IndexDelimitedLocusShard next() {
        FilePointer nextFilePointer = filePointerIterator.next();
        Map<SAMFileReader2,List<Chunk>> chunksBounding = dataSource.getFilePointersBounding(nextFilePointer.bin);
        return new IndexDelimitedLocusShard(nextFilePointer.locations,chunksBounding,Shard.ShardType.LOCUS_INTERVAL);
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