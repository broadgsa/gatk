/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package htsjdk.samtools;

import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

/**
 * A temporary solution to work around Java access rights issues:
 * override BAMFileSpan and make it public.
 * TODO: Eliminate once we determine the final fate of the BAM index reading code.
 */
public class GATKBAMFileSpan extends BAMFileSpan {
    /**
     * Create a new empty list of chunks.
     */
    public GATKBAMFileSpan() {
        super();
    }

    /**
     * Create a new GATKBAMFileSpan from an existing BAMFileSpan.
     * @param sourceFileSpan
     */
    public GATKBAMFileSpan(SAMFileSpan sourceFileSpan) {
        if(!(sourceFileSpan instanceof BAMFileSpan))
            throw new SAMException("Unable to create GATKBAMFileSpan from a SAMFileSpan. Please submit a BAMFileSpan instead");
        BAMFileSpan sourceBAMFileSpan = (BAMFileSpan)sourceFileSpan;
        for(Chunk chunk: sourceBAMFileSpan.getChunks())
            add(chunk instanceof GATKChunk ? chunk : new GATKChunk(chunk));
    }

    /**
     * Convenience constructor to construct a BAM file span from
     * a single chunk.
     * @param chunk Chunk to use as the sole region in this span.
     */
    public GATKBAMFileSpan(final Chunk chunk) {
        super(chunk);
    }

    /**
     * Create a new chunk list from the given list of chunks.
     * @param chunks Constituent chunks.
     */
    public GATKBAMFileSpan(final GATKChunk[] chunks) {
        super(Arrays.<Chunk>asList(chunks));
    }

    @Override
    public boolean equals(final Object other) {
        if(!(other instanceof BAMFileSpan))
            return false;

        List<Chunk> theseChunks = getChunks();
        List<Chunk> otherChunks = ((BAMFileSpan)other).getChunks();

        if(theseChunks.size() != otherChunks.size())
            return false;
        for(int i = 0; i < theseChunks.size(); i++) {
            if(!theseChunks.get(i).equals(otherChunks.get(i)))
                return false;
        }

        return true;
    }

    /**
     * Gets the constituent chunks stored in this span.
     * @return An unmodifiable list of chunks.
     */
    public List<GATKChunk> getGATKChunks() {
        List<GATKChunk> gatkChunks = new ArrayList<GATKChunk>();
        for(Chunk chunk: getChunks())
            gatkChunks.add(new GATKChunk(chunk));
        return gatkChunks;
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        for(GATKChunk chunk: getGATKChunks())
            builder.append(String.format("%s;",chunk));
        return builder.toString();
    }

    /**
     * Returns an approximation of the number of uncompressed bytes in this
     * file span.
     * @return Approximation of uncompressed bytes in filespan.
     */
    public long size() {
        long size = 0L;
        for(GATKChunk chunk: getGATKChunks())
            size += chunk.size();
        return size;
    }

    /**
     * Get a GATKChunk representing the "extent" of this file span, from the start of the first
     * chunk to the end of the last chunk.The chunks list must be sorted in order to use this method.
     *
     * @return a GATKChunk representing the extent of this file span, or a GATKChunk representing
     *         a span of size 0 if there are no chunks
     */
    public GATKChunk getExtent() {
        validateSorted();   // TODO: defensive measure: may be unnecessary

        List<Chunk> chunks = getChunks();
        if ( chunks.isEmpty() ) {
            return new GATKChunk(0L, 0L);
        }

        return new GATKChunk(chunks.get(0).getChunkStart(), chunks.get(chunks.size() - 1).getChunkEnd());
    }

    /**
     * Validates the list of chunks to ensure that they appear in sorted order.
     */
    private void validateSorted() {
        List<Chunk> chunks = getChunks();
        for ( int i = 1; i < chunks.size(); i++ ) {
            if ( chunks.get(i).getChunkStart() < chunks.get(i-1).getChunkEnd() ) {
                throw new ReviewedGATKException(String.format("Chunk list is unsorted; chunk %s is before chunk %s", chunks.get(i-1), chunks.get(i)));

            }
        }
    }

    /**
     * Computes the union of two FileSpans.
     * @param other FileSpan to union with this one.
     * @return A file span that's been unioned.
     */
    public GATKBAMFileSpan union(final GATKBAMFileSpan other) {
        // No data?  Return an empty file span.
        if(getGATKChunks().size() == 0 && other.getGATKChunks().size() == 0)
            return new GATKBAMFileSpan();

        LinkedList<GATKChunk> unmergedUnion = new LinkedList<GATKChunk>();
        unmergedUnion.addAll(getGATKChunks());
        unmergedUnion.addAll(other.getGATKChunks());
        Collections.sort(unmergedUnion);

        List<GATKChunk> mergedUnion = new ArrayList<GATKChunk>();
        GATKChunk currentChunk = unmergedUnion.remove();
        while(!unmergedUnion.isEmpty()) {

            // While the current chunk can be merged with the next chunk:
            while( ! unmergedUnion.isEmpty() &&
                   (currentChunk.overlaps(unmergedUnion.peek()) || currentChunk.isAdjacentTo(unmergedUnion.peek())) ) {

                // Merge the current chunk with the next chunk:
                GATKChunk nextChunk = unmergedUnion.remove();
                currentChunk = currentChunk.merge(nextChunk);
            }
            // Add the accumulated range.
            mergedUnion.add(currentChunk);
            currentChunk = !unmergedUnion.isEmpty() ? unmergedUnion.remove() : null;
        }

        // At end of the loop above, the last chunk will be contained in currentChunk and will not yet have been added.  Add it.
        if(currentChunk !=null)
            mergedUnion.add(currentChunk);

        return new GATKBAMFileSpan(mergedUnion.toArray(new GATKChunk[mergedUnion.size()]));
    }

    /**
     * Intersects two BAM file spans.
     * @param other File span to intersect with this one.
     * @return The intersected BAM file span.
     */
    public GATKBAMFileSpan intersection(final GATKBAMFileSpan other) {
        Iterator<GATKChunk> thisIterator = getGATKChunks().iterator();
        Iterator<GATKChunk> otherIterator = other.getGATKChunks().iterator();

        if(!thisIterator.hasNext() || !otherIterator.hasNext())
            return new GATKBAMFileSpan();

        GATKChunk thisChunk = thisIterator.next();
        GATKChunk otherChunk = otherIterator.next();

        List<GATKChunk> intersected = new ArrayList<GATKChunk>();

        while(thisChunk != null && otherChunk != null) {
            // If this iterator is before other, skip this ahead.
            if(thisChunk.getChunkEnd() <= otherChunk.getChunkStart()) {
                thisChunk = thisIterator.hasNext() ? thisIterator.next() : null;
                continue;
            }

            // If other iterator is before this, skip other ahead.
            if(thisChunk.getChunkStart() >= otherChunk.getChunkEnd()) {
                otherChunk = otherIterator.hasNext() ? otherIterator.next() : null;
                continue;
            }

            // If these two chunks overlap, pull out intersection of data and truncated current chunks to point after
            // the intersection (or next chunk if no such overlap exists).
            if(thisChunk.overlaps(otherChunk)) {
                // Determine the chunk constraints
                GATKChunk firstChunk = thisChunk.getChunkStart() < otherChunk.getChunkStart() ? thisChunk : otherChunk;
                GATKChunk secondChunk = thisChunk==firstChunk ? otherChunk : thisChunk;
                GATKChunk intersectedChunk = new GATKChunk(secondChunk.getChunkStart(),Math.min(firstChunk.getChunkEnd(),secondChunk.getChunkEnd()));
                intersected.add(intersectedChunk);

                if(thisChunk.getChunkEnd() > intersectedChunk.getChunkEnd())
                    thisChunk = new GATKChunk(intersectedChunk.getChunkEnd(),thisChunk.getChunkEnd());
                else
                    thisChunk = thisIterator.hasNext() ? thisIterator.next() : null;

                if(otherChunk.getChunkEnd() > intersectedChunk.getChunkEnd())
                    otherChunk = new GATKChunk(intersectedChunk.getChunkEnd(),otherChunk.getChunkEnd());
                else
                    otherChunk = otherIterator.hasNext() ? otherIterator.next() : null;
            }

        }

        return new GATKBAMFileSpan(intersected.toArray(new GATKChunk[intersected.size()]));
    }

    /**
     * Substracts other file span from this file span.
     * @param other File span to strike out.
     * @return This file span minuse the other file span.
     */

    public GATKBAMFileSpan minus(final GATKBAMFileSpan other) {
        Iterator<GATKChunk> thisIterator = getGATKChunks().iterator();
        Iterator<GATKChunk> otherIterator = other.getGATKChunks().iterator();

        if(!thisIterator.hasNext() || !otherIterator.hasNext())
            return this;

        GATKChunk thisChunk = thisIterator.next();
        GATKChunk otherChunk = otherIterator.next();

        List<GATKChunk> subtracted = new ArrayList<GATKChunk>();

        while(thisChunk != null && otherChunk != null) {
            // If this iterator is before the other, add this to the subtracted list and forge ahead.
            if(thisChunk.getChunkEnd() <= otherChunk.getChunkStart()) {
                subtracted.add(thisChunk);
                thisChunk = thisIterator.hasNext() ? thisIterator.next() : null;
                continue;
            }

            // If other iterator is before this, skip other ahead.
            if(thisChunk.getChunkStart() >= otherChunk.getChunkEnd()) {
                otherChunk = otherIterator.hasNext() ? otherIterator.next() : null;
                continue;
            }

            // If these two chunks overlap, pull out intersection of data and truncated current chunks to point after
            // the intersection (or next chunk if no such overlap exists).
            if(thisChunk.overlaps(otherChunk)) {
                // Add in any sort of prefix that this chunk might have over the other.
                if(thisChunk.getChunkStart() < otherChunk.getChunkStart())
                    subtracted.add(new GATKChunk(thisChunk.getChunkStart(),otherChunk.getChunkStart()));

                if(thisChunk.getChunkEnd() > otherChunk.getChunkEnd())
                    thisChunk = new GATKChunk(otherChunk.getChunkEnd(),thisChunk.getChunkEnd());
                else
                    thisChunk = thisIterator.hasNext() ? thisIterator.next() : null;
            }
        }

        // Finish up any remaining contents of this that didn't make it into the subtracted array.
        if(thisChunk != null)
            subtracted.add(thisChunk);
        while(thisIterator.hasNext())
            subtracted.add(thisIterator.next());

        return new GATKBAMFileSpan(subtracted.toArray(new GATKChunk[subtracted.size()]));
    }
}
