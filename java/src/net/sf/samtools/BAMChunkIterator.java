package net.sf.samtools;

import java.util.*;

/**
 * Adapts a list of BAM blocks to BAM chunks.
 *
 * @author mhanna
 * @version 0.1
 */
public class BAMChunkIterator implements Iterator<Chunk> {
    /**
     * The wrapped chunk iterator.
     */
    private final BAMBlockIterator blockIterator;

    /**
     * List of prefetched segments.  If prefetchedSegments is
     * empty, this should mean that there's nothing left in the queue.
     */
    private Queue<BlockSegment> prefetchedSegments;

    /**
     * List of chunks to filter out of the final chunk list.
     */
    private final Queue<Chunk> filters;

    /**
     * Creates a new chunk iterator wrapping the given block, optionally filtering
     * out a list of chunks from those blocks.
     * TODO: This class is too smart.  Wrap filtering functionality in another class.
     * @param blockIterator block iterator to wrap.
     * @param filters List of chunks to filter out of the given input stream.
     */
    public BAMChunkIterator(BAMBlockIterator blockIterator, List<Chunk> filters) {
        this.blockIterator = blockIterator;
        this.prefetchedSegments = new LinkedList<BlockSegment>();
        this.filters = new PriorityQueue<Chunk>(filters);
        seedNextSegments();
    }

    /**
     * Is a next chunk available.
     * @return true if a next chunk is available; false otherwise.
     */
    public boolean hasNext() {
        return !prefetchedSegments.isEmpty();    
    }

    /**
     * Returns the next chunk in the list.
     * @return Next block, adapted to a chunk.
     */
    public Chunk next() {
        if(prefetchedSegments.isEmpty())
            throw new NoSuchElementException("BAMChunkIterator is exhausted.");
        // TODO: Maybe tie together adjacent chunks together?
        Chunk nextChunk = prefetchedSegments.remove().toChunk();
        seedNextSegments();
        return nextChunk;
    }

    /**
     * @throws UnsupportedOperationException always.
     */
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove from a BAM chunk iterator");
    }

    /**
     * Seed the next value or set of values for the iterator, if necessary.
     */
    private void seedNextSegments() {
        while(prefetchedSegments.isEmpty() && blockIterator.hasNext()) {
            // Add the new segments to the prefetch...
            prefetchedSegments.add(new BlockSegment(blockIterator.next()));

            // And iteratively subtract the filtering chunks.
            for(Chunk filter: filters) {
                Queue<BlockSegment> filteredBlockSegments = new LinkedList<BlockSegment>();
                for(BlockSegment blockSegment: prefetchedSegments)
                    filteredBlockSegments.addAll(blockSegment.minus(filter));
                prefetchedSegments = filteredBlockSegments;
            }
        }
    }

}
