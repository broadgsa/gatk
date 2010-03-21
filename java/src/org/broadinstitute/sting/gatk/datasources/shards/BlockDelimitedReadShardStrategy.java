package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.*;
import net.sf.picard.filter.SamRecordFilter;

import java.util.*;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

/**
 * A read shard strategy that delimits based on the number of
 * blocks in the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class BlockDelimitedReadShardStrategy extends ReadShardStrategy {
    /**
     * What is the maximum number of reads which should go into a read shard.
     */
    protected static final int MAX_READS = 10000;

    /**
     * The data source used to shard.
     */
    protected final BlockDrivenSAMDataSource dataSource;

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
    private Map<SAMReaderID,Chunk> position;

    /**
     * Create a new read shard strategy, loading read shards from the given BAM file.
     * @param dataSource Data source from which to load shards.
     * @param locations intervals to use for sharding.
     */
    public BlockDelimitedReadShardStrategy(SAMDataSource dataSource, GenomeLocSortedSet locations) {
        if(!(dataSource instanceof BlockDrivenSAMDataSource))
            throw new StingException("Block-delimited read shard strategy cannot work without a block-driven data source.");

        this.dataSource = (BlockDrivenSAMDataSource)dataSource;
        this.position = this.dataSource.getCurrentPosition();
        if(locations != null)
            filePointers.addAll(IntervalSharder.shardIntervals(this.dataSource,locations.toList()));

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
     * @throws NoSuchElementException if no such shard is available.
     */
    public Shard next() {
        if(!hasNext())
            throw new NoSuchElementException("No next read shard available");
        Shard currentShard = nextShard;
        advance();
        return currentShard;
    }

    public void advance() {
        Map<SAMReaderID,List<Chunk>> shardPosition = new HashMap<SAMReaderID,List<Chunk>>();
        nextShard = null;
        SamRecordFilter filter = null;

        if(!filePointers.isEmpty()) {
            Map<SAMReaderID,List<Chunk>> selectedReaders = new HashMap<SAMReaderID,List<Chunk>>();
            while(selectedReaders.size() == 0 && currentFilePointer != null) {
                shardPosition = currentFilePointer.chunks;
                for(SAMReaderID id: shardPosition.keySet()) {
                    List<Chunk> chunks = shardPosition.get(id);
                    List<Chunk> selectedChunks = new ArrayList<Chunk>();
                    Chunk filePosition = position.get(id);
                    for(Chunk chunk: chunks)
                        if(filePosition.getChunkStart() <= chunk.getChunkStart())
                            selectedChunks.add(chunk);
                        else if(filePosition.getChunkStart() > chunk.getChunkStart() && filePosition.getChunkStart() < chunk.getChunkEnd()) {
                            selectedChunks.add(new Chunk(filePosition.getChunkStart(),chunk.getChunkEnd()));
                    }
                    if(selectedChunks.size() > 0)
                        selectedReaders.put(id,selectedChunks);
                }

                if(selectedReaders.size() > 0) {
                    filter = new ReadOverlapFilter(currentFilePointer.locations);
                    BAMFormatAwareShard shard = new BlockDelimitedReadShard(dataSource.getReadsInfo(),selectedReaders,filter,Shard.ShardType.READ);
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
            // TODO: This level of processing should not be necessary.
            shardPosition = new HashMap<SAMReaderID,List<Chunk>>();
            for(Map.Entry<SAMReaderID,Chunk> entry: position.entrySet())
                shardPosition.put(entry.getKey(),Collections.singletonList(entry.getValue()));
            filter = null;

            BAMFormatAwareShard shard = new BlockDelimitedReadShard(dataSource.getReadsInfo(),shardPosition,filter,Shard.ShardType.READ);
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
