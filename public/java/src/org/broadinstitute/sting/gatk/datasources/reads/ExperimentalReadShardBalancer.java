/*
 * Copyright (c) 2012, The Broad Institute
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

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMFileSpan;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;

import java.util.*;

/**
 * Convert from an unbalanced iterator over FilePointers to a balanced iterator over Shards
 *
 * @author David Roazen
 */
public class ExperimentalReadShardBalancer extends ShardBalancer {

    private static Logger logger = Logger.getLogger(ExperimentalReadShardBalancer.class);

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
            private FilePointer currentFilePointer = null;

            /**
             * Iterator over the reads from the current file pointer. The same iterator will be
             * used to fill all shards associated with a given file pointer
             */
            private PeekableIterator<SAMRecord> currentFilePointerReadsIterator = null;

            {
                if ( filePointers.hasNext() )
                    currentFilePointer = filePointers.next();
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
                while ( nextShard == null && currentFilePointer != null ) {

                    // If we've exhausted the current file pointer of reads, move to the next file pointer (if there is one):
                    if ( currentFilePointerReadsIterator != null && ! currentFilePointerReadsIterator.hasNext() ) {

                        // Close the old, exhausted chain of iterators to release resources
                        currentFilePointerReadsIterator.close();

                        do {
                            advanceFilePointer();
                        } while ( currentFilePointer != null && isEmpty(currentFilePointer.fileSpans) ); // skip empty file pointers

                        // We'll need to create a fresh iterator for this file pointer when we create the first
                        // shard for it below.
                        currentFilePointerReadsIterator = null;
                    }

                    // At this point if currentFilePointer is non-null we know it is also non-empty. Our
                    // currentFilePointerReadsIterator may be null or non-null depending on whether or not
                    // this is our first shard for this file pointer.
                    if ( currentFilePointer != null ) {
                        Shard shard = new ReadShard(parser,readsDataSource,currentFilePointer.fileSpans,currentFilePointer.locations,currentFilePointer.isRegionUnmapped);

                        // Create a new reads iterator only when we've just advanced to a new file pointer. It's
                        // essential that the iterators persist across all shards that share the same file pointer
                        // to allow the downsampling to work properly.
                        if ( currentFilePointerReadsIterator == null ) {
                            currentFilePointerReadsIterator = new PeekableIterator<SAMRecord>(readsDataSource.getIterator(shard));
                        }

                        if ( currentFilePointerReadsIterator.hasNext() ) {
                            shard.fill(currentFilePointerReadsIterator);
                            nextShard = shard;
                        }
                    }
                }
            }

            private void advanceFilePointer() {
                FilePointer previousFilePointer = currentFilePointer;
                currentFilePointer = filePointers.hasNext() ? filePointers.next() : null;

                // TODO: This is a purely defensive measure to guard against the possibility of overlap
                // TODO: between FilePointers. When overlap is detected, remove the overlapping regions from
                // TODO: the newly-current FilePointer.
                // TODO: If we later discover that overlap is theoretically impossible, this step becomes
                // TODO: unnecessary and should be removed.
                if ( currentFilePointer != null && previousFilePointer != null &&
                     previousFilePointer.hasFileSpansOverlappingWith(currentFilePointer) ) {

                    logger.debug(String.format("%s: found consecutive overlapping FilePointers [%s] and [%s]", getClass().getSimpleName(), previousFilePointer, currentFilePointer));

                    Map<SAMReaderID, SAMFileSpan> previousFileSpans = previousFilePointer.getFileSpans();
                    Map<SAMReaderID, SAMFileSpan> trimmedFileSpans = new HashMap<SAMReaderID, SAMFileSpan>(currentFilePointer.getFileSpans().size());

                    for ( Map.Entry<SAMReaderID, SAMFileSpan> fileSpanEntry : currentFilePointer.getFileSpans().entrySet() ) {
                        // find the corresponding file span from the previous FilePointer
                        SAMFileSpan previousFileSpan = previousFileSpans.get(fileSpanEntry.getKey());

                        if ( previousFileSpan == null ) {
                            // no match, so no trimming required
                            trimmedFileSpans.put(fileSpanEntry.getKey(), fileSpanEntry.getValue());
                        }
                        else {
                            // match, so remove any overlapping regions (regions before the start of the
                            // region immediately following the previous file span)
                            SAMFileSpan trimmedSpan = fileSpanEntry.getValue().removeContentsBefore(previousFileSpan.getContentsFollowing());
                            trimmedFileSpans.put(fileSpanEntry.getKey(), trimmedSpan);
                        }
                    }

                    // Replace the current file pointer with its trimmed equivalent
                    currentFilePointer = new FilePointer(trimmedFileSpans, currentFilePointer.locations);
                }
            }

            /**
             * Detects whether the list of file spans contain any read data.
             * @param selectedSpans Mapping of readers to file spans.
             * @return True if file spans are completely empty; false otherwise.
             */
            private boolean isEmpty(Map<SAMReaderID,SAMFileSpan> selectedSpans) {
                for(SAMFileSpan fileSpan: selectedSpans.values()) {
                    if(!fileSpan.isEmpty())
                        return false;
                }
                return true;
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from shard balancing iterator");
            }
        };
    }

}
