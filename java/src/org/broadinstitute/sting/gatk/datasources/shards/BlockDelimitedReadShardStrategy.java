package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.samtools.*;

import java.util.*;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;

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

    /** our storage of the genomic locations they'd like to shard over */
    private final List<FilePointer> filePointers = new ArrayList<FilePointer>();

    /**
     * Position of the last shard in the file.
     */
    private Map<SAMFileReader2,Chunk> position;

    /**
     * True if the set of SAM readers is at the end of the stream.  False otherwise.
     */
    private boolean atEndOfStream = false;

    /**
     * Create a new read shard strategy, loading read shards from the given BAM file.
     * @param dataSource Data source from which to load shards.
     */
    public BlockDelimitedReadShardStrategy(SAMDataSource dataSource, GenomeLocSortedSet locations) {
        if(!(dataSource instanceof BlockDrivenSAMDataSource))
            throw new StingException("Block-delimited read shard strategy cannot work without a block-driven data source.");

        this.dataSource = (BlockDrivenSAMDataSource)dataSource;
        this.position = this.dataSource.getCurrentPosition();
        if(locations != null)
            filePointers.addAll(IntervalSharder.shardIntervals(this.dataSource,locations.toList(),this.dataSource.getNumIndexLevels()-1));
    }

    /**
     * do we have another read shard?
     * @return True if any more data is available.  False otherwise.
     */
    public boolean hasNext() {
        return !atEndOfStream;
    }

    /**
     * Retrieves the next shard, if available.
     * @return The next shard, if available.
     * @throws NoSuchElementException if no such shard is available.
     */
    public Shard next() {
        if(!hasNext())
            throw new NoSuchElementException("No such element available: SAM reader has arrived at last shard.");

        Map<SAMFileReader2,List<Chunk>> shardPosition = null;
        if(!filePointers.isEmpty()) {
            boolean foundData = false;
            for(FilePointer filePointer: filePointers) {
                shardPosition = dataSource.getFilePointersBounding(filePointer.bin);
                for(SAMFileReader2 reader: shardPosition.keySet()) {
                    List<Chunk> chunks = shardPosition.get(reader);
                    Chunk filePosition = position.get(reader);
                    for(Chunk chunk: chunks) {
                        if(filePosition.getChunkStart() > chunk.getChunkEnd())
                            chunks.remove(chunk);
                        else {
                            if(filePosition.getChunkStart() > chunk.getChunkStart())
                                chunk.setChunkStart(filePosition.getChunkStart());
                            foundData = true;
                        }
                    }
                }
                if(foundData)
                    break;
            }
        }
        else {
            // TODO: This level of processing should not be necessary.
            shardPosition = new HashMap<SAMFileReader2,List<Chunk>>();
            for(Map.Entry<SAMFileReader2,Chunk> entry: position.entrySet())
                shardPosition.put(entry.getKey(),Collections.singletonList(entry.getValue()));
        }

        BAMFormatAwareShard shard = new BlockDelimitedReadShard(dataSource.getReadsInfo(),shardPosition,Shard.ShardType.READ);
        atEndOfStream = dataSource.fillShard(shard);

        this.position = dataSource.getCurrentPosition();

        return shard;        
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
