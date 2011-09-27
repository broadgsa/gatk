/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.samtools.SAMFileSpan;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.*;

/**
 * The sharding strategy for reads using a simple counting mechanism.  Each read shard
 * has a specific number of reads (default to 10K) which is configured in the constructor.
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 */
public class ReadShardStrategy implements ShardStrategy {
    /**
     * What is the maximum number of reads which should go into a read shard.
     */
    protected static final int MAX_READS = 10000;

    /**
     * The data source used to shard.
     */
    private final SAMDataSource dataSource;

    /**
     * The intervals to be processed.
     */
    private final GenomeLocSortedSet locations;

    /**
     * The cached shard to be returned next.  Prefetched in the peekable iterator style.
     */
    private Shard nextShard = null;

    /** our storage of the genomic locations they'd like to shard over */
    private final List<FilePointer> filePointers = new ArrayList<FilePointer>();

    /**
     * Iterator over the list of file pointers.
     */
    private final Iterator<FilePointer> filePointerIterator;

    /**
     * The file pointer currently being processed.
     */
    private FilePointer currentFilePointer;

    /**
     * Ending position of the last shard in the file.
     */
    private Map<SAMReaderID,SAMFileSpan> position;

    /**
     * An indicator whether the strategy has sharded into the unmapped region.
     */
    private boolean isIntoUnmappedRegion = false;

    private final GenomeLocParser parser;

    /**
     * Create a new read shard strategy, loading read shards from the given BAM file.
     * @param dataSource Data source from which to load shards.
     * @param locations intervals to use for sharding.
     */
    public ReadShardStrategy(GenomeLocParser parser, SAMDataSource dataSource, GenomeLocSortedSet locations) {
        this.dataSource = dataSource;
        this.parser = parser;
        this.position = this.dataSource.getCurrentPosition();
        this.locations = locations;

        if(locations != null)
            filePointerIterator = dataSource.isLowMemoryShardingEnabled() ? new LowMemoryIntervalSharder(this.dataSource,locations) : IntervalSharder.shardIntervals(this.dataSource,locations);
        else
            filePointerIterator = filePointers.iterator();

        if(filePointerIterator.hasNext())
            currentFilePointer = filePointerIterator.next();

        advance();
    }

    /**
     * do we have another read shard?
     * @return True if any more data is available.  False otherwise.
     */
    public boolean hasNext() {
        return nextShard != null;
    }

    /**
     * Retrieves the next shard, if available.
     * @return The next shard, if available.
     * @throws java.util.NoSuchElementException if no such shard is available.
     */
    public Shard next() {
        if(!hasNext())
            throw new NoSuchElementException("No next read shard available");
        Shard currentShard = nextShard;
        advance();
        return currentShard;
    }

    public void advance() {
        Map<SAMReaderID,SAMFileSpan> shardPosition = new HashMap<SAMReaderID,SAMFileSpan>();
        nextShard = null;

        if(locations != null) {
            Map<SAMReaderID,SAMFileSpan> selectedReaders = new HashMap<SAMReaderID,SAMFileSpan>();
            while(selectedReaders.size() == 0 && currentFilePointer != null) {
                shardPosition = currentFilePointer.fileSpans;

                for(SAMReaderID id: shardPosition.keySet()) {
                    SAMFileSpan fileSpan = shardPosition.get(id).removeContentsBefore(position.get(id));
                    if(!fileSpan.isEmpty())
                        selectedReaders.put(id,fileSpan);
                }

                if(selectedReaders.size() > 0) {
                    Shard shard = new ReadShard(parser, dataSource,selectedReaders,currentFilePointer.locations,currentFilePointer.isRegionUnmapped);
                    dataSource.fillShard(shard);

                    if(!shard.isBufferEmpty()) {
                        nextShard = shard;
                        break;
                    }
                }

                selectedReaders.clear();
                currentFilePointer = filePointerIterator.hasNext() ? filePointerIterator.next() : null;
            }
        }
        else {
            // todo -- this nulling of intervals is a bit annoying since readwalkers without
            // todo -- any -L values need to be special cased throughout the code.
            Shard shard = new ReadShard(parser,dataSource,position,null,false);
            dataSource.fillShard(shard);
            nextShard = !shard.isBufferEmpty() ? shard : null;
        }

        this.position = dataSource.getCurrentPosition();
    }

    /**
     * @throws UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Remove not supported");
    }

    /**
     * Convenience method for using ShardStrategy in an foreach loop.
     * @return A iterator over shards.
     */
    public Iterator<Shard> iterator() {
        return this;
    }
}
