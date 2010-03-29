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
    private final BlockDrivenSAMDataSource dataSource;

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
    private Map<SAMReaderID,BAMFileSpan> position;

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
        this.locations = locations;

        if(locations != null)
            filePointerIterator = IntervalSharder.shardIntervals(this.dataSource,locations.toList());
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
        Map<SAMReaderID,BAMFileSpan> shardPosition = new HashMap<SAMReaderID,BAMFileSpan>();
        nextShard = null;
        SamRecordFilter filter = null;

        if(locations != null) {
            Map<SAMReaderID,BAMFileSpan> selectedReaders = new HashMap<SAMReaderID,BAMFileSpan>();
            while(selectedReaders.size() == 0 && currentFilePointer != null) {
                shardPosition = currentFilePointer.fileSpans;
                for(SAMReaderID id: shardPosition.keySet()) {
                    BAMFileSpan fileSpans = shardPosition.get(id).removeBefore(position.get(id));
                    if(!fileSpans.isEmpty())
                        selectedReaders.put(id,fileSpans);
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
            BAMFormatAwareShard shard = new BlockDelimitedReadShard(dataSource.getReadsInfo(),position,filter,Shard.ShardType.READ);
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
