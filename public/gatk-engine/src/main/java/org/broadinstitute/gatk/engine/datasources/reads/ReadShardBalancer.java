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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * Convert from an unbalanced iterator over FilePointers to a balanced iterator over Shards.
 *
 * When processing FilePointers, our strategy is to aggregate all FilePointers for each contig
 * together into one monolithic FilePointer, create one persistent set of read iterators over
 * that monolithic FilePointer, and repeatedly use that persistent set of read iterators to
 * fill read shards with reads.
 *
 * This strategy has several important advantages:
 *
 * 1. We avoid issues with file span overlap. FilePointers that are more granular than a whole
 *    contig will have regions that overlap with other FilePointers on the same contig, due
 *    to the limited granularity of BAM index data. By creating only one FilePointer per contig,
 *    we avoid having to track how much of each file region we've visited (as we did in the
 *    former implementation), we avoid expensive non-sequential access patterns in the files,
 *    and we avoid having to repeatedly re-create our iterator chain for every small region
 *    of interest.
 *
 * 2. We avoid boundary issues with the engine-level downsampling. Since we create a single
 *    persistent set of read iterators (which include the downsampling iterator(s)) per contig,
 *    the downsampling process is never interrupted by FilePointer or Shard boundaries, and never
 *    loses crucial state information while downsampling within a contig.
 *
 * TODO: There is also at least one important disadvantage:
 *
 * 1. We load more BAM index data into memory at once, and this work is done upfront before processing
 *    the next contig, creating a delay before traversal of each contig. This delay may be
 *    compensated for by the gains listed in #1 above, and we may be no worse off overall in
 *    terms of total runtime, but we need to verify this empirically.
 *
 * @author David Roazen
 */
public class ReadShardBalancer extends ShardBalancer {

    private static Logger logger = Logger.getLogger(ReadShardBalancer.class);

    /**
     * Convert iterators of file pointers into balanced iterators of shards.
     * @return An iterator over balanced shards.
     */
    public Iterator<Shard> iterator() {
        return new Iterator<Shard>() {
            /**
             * The cached shard to be returned next.  Prefetched in the peekable iterator style.
             */
            private Shard nextShard = null;

            /**
             * The file pointer currently being processed.
             */
            private FilePointer currentContigFilePointer = null;

            /**
             * Iterator over the reads from the current contig's file pointer. The same iterator will be
             * used to fill all shards associated with a given file pointer
             */
            private PeekableIterator<SAMRecord> currentContigReadsIterator = null;

            /**
             * How many FilePointers have we pulled from the filePointers iterator?
             */
            private int totalFilePointersConsumed = 0;

            /**
             * Have we encountered a monolithic FilePointer?
             */
            private boolean encounteredMonolithicFilePointer = false;


            {
                createNextContigFilePointer();
                advance();
            }

            public boolean hasNext() {
                return nextShard != null;
            }

            public Shard next() {
                if ( ! hasNext() )
                    throw new NoSuchElementException("No next read shard available");
                Shard currentShard = nextShard;
                advance();
                return currentShard;
            }

            private void advance() {
                nextShard = null;

                // May need multiple iterations to fill the next shard if all reads in current file spans get filtered/downsampled away
                while ( nextShard == null && currentContigFilePointer != null ) {

                    // If we've exhausted the current file pointer of reads, move to the next file pointer (if there is one):
                    if ( currentContigReadsIterator != null && ! currentContigReadsIterator.hasNext() ) {

                        // Close the old, exhausted chain of iterators to release resources
                        currentContigReadsIterator.close();

                        // Advance to the FilePointer for the next contig
                        createNextContigFilePointer();

                        // We'll need to create a fresh iterator for this file pointer when we create the first
                        // shard for it below.
                        currentContigReadsIterator = null;
                    }

                    // At this point our currentContigReadsIterator may be null or non-null depending on whether or not
                    // this is our first shard for this file pointer.
                    if ( currentContigFilePointer != null ) {
                        Shard shard = new ReadShard(parser,readsDataSource, currentContigFilePointer.fileSpans, currentContigFilePointer.locations, currentContigFilePointer.isRegionUnmapped);

                        // Create a new reads iterator only when we've just advanced to the file pointer for the next
                        // contig. It's essential that the iterators persist across all shards that share the same contig
                        // to allow the downsampling to work properly.
                        if ( currentContigReadsIterator == null ) {
                            currentContigReadsIterator = new PeekableIterator<SAMRecord>(readsDataSource.getIterator(shard));
                        }

                        if ( currentContigReadsIterator.hasNext() ) {
                            shard.fill(currentContigReadsIterator);
                            nextShard = shard;
                        }
                    }
                }
            }

            /**
             * Aggregate all FilePointers for the next contig together into one monolithic FilePointer
             * to avoid boundary issues with visiting the same file regions more than once (since more
             * granular FilePointers will have regions that overlap with other nearby FilePointers due
             * to the nature of BAM indices).
             *
             * By creating one persistent set of iterators per contig we also avoid boundary artifacts
             * in the engine-level downsampling.
             *
             * TODO: This FilePointer aggregation should ideally be done at the BAMSchedule level for
             * TODO: read traversals, as there's little point in the BAMSchedule emitting extremely
             * TODO: granular FilePointers if we're just going to union them. The BAMSchedule should
             * TODO: emit one FilePointer per contig for read traversals (but, crucially, NOT for
             * TODO: locus traversals).
             */
            private void createNextContigFilePointer() {
                currentContigFilePointer = null;
                List<FilePointer> nextContigFilePointers = new ArrayList<FilePointer>();

                if ( filePointers.hasNext() ) {
                    logger.info("Loading BAM index data");
                }

                while ( filePointers.hasNext() ) {

                    // Make sure that if we see a monolithic FilePointer (representing all regions in all files) that
                    // it is the ONLY FilePointer we ever encounter
                    if ( encounteredMonolithicFilePointer ) {
                        throw new ReviewedGATKException("Bug: encountered additional FilePointers after encountering a monolithic FilePointer");
                    }
                    if ( filePointers.peek().isMonolithic() ) {
                        if ( totalFilePointersConsumed > 0 ) {
                            throw new ReviewedGATKException("Bug: encountered additional FilePointers before encountering a monolithic FilePointer");
                        }
                        encounteredMonolithicFilePointer = true;
                        logger.debug(String.format("Encountered monolithic FilePointer: %s", filePointers.peek()));
                    }

                    // If this is the first FP we've seen, or we're dealing with mapped regions and the next FP is on the
                    // same contig as previous FPs, or all our FPs are unmapped, add the next FP to the list of FPs to merge
                    if ( nextContigFilePointers.isEmpty() ||
                             (! nextContigFilePointers.get(0).isRegionUnmapped && ! filePointers.peek().isRegionUnmapped &&
                             nextContigFilePointers.get(0).getContigIndex() == filePointers.peek().getContigIndex()) ||
                                 (nextContigFilePointers.get(0).isRegionUnmapped && filePointers.peek().isRegionUnmapped) ) {

                        nextContigFilePointers.add(filePointers.next());
                        totalFilePointersConsumed++;
                    }
                    else {
                        break; // next FilePointer is on a different contig or has different mapped/unmapped status,
                               // save it for next time
                    }
                }

                if ( ! nextContigFilePointers.isEmpty() ) {
                    currentContigFilePointer = FilePointer.union(nextContigFilePointers, parser);
                }

                if ( currentContigFilePointer != null ) {
                    logger.info("Done loading BAM index data");
                    logger.debug(String.format("Next FilePointer: %s", currentContigFilePointer));
                }
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from shard balancing iterator");
            }
        };
    }

}
