/*
 * Copyright (c) 2011, The Broad Institute
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

import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.SAMFileSpan;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

/**
 * Divide up large file pointers containing reads into more manageable subcomponents.
 */
public class ReadShardBalancer extends ShardBalancer {
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
            private FilePointer currentFilePointer;

            /**
             * Ending position of the last shard in the file.
             */
            private Map<SAMReaderID,GATKBAMFileSpan> position = readsDataSource.getCurrentPosition();

            {
                if(filePointers.hasNext())
                    currentFilePointer = filePointers.next();
                advance();
            }

            public boolean hasNext() {
                return nextShard != null;
            }

            public Shard next() {
                if(!hasNext())
                    throw new NoSuchElementException("No next read shard available");
                Shard currentShard = nextShard;
                advance();
                return currentShard;
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from shard balancing iterator");
            }

            private void advance() {
                Map<SAMReaderID,SAMFileSpan> shardPosition;
                nextShard = null;

                Map<SAMReaderID,SAMFileSpan> selectedReaders = new HashMap<SAMReaderID,SAMFileSpan>();
                while(selectedReaders.size() == 0 && currentFilePointer != null) {
                    shardPosition = currentFilePointer.fileSpans;

                    for(SAMReaderID id: shardPosition.keySet()) {
                        SAMFileSpan fileSpan = new GATKBAMFileSpan(shardPosition.get(id).removeContentsBefore(position.get(id)));
                        if(!fileSpan.isEmpty())
                            selectedReaders.put(id,fileSpan);
                    }

                    if(selectedReaders.size() > 0) {
                        Shard shard = new ReadShard(parser,readsDataSource,selectedReaders,currentFilePointer.locations,currentFilePointer.isRegionUnmapped);
                        readsDataSource.fillShard(shard);

                        if(!shard.isBufferEmpty()) {
                            nextShard = shard;
                            break;
                        }
                    }

                    selectedReaders.clear();
                    currentFilePointer = filePointers.hasNext() ? filePointers.next() : null;
                }

                position = readsDataSource.getCurrentPosition();
            }
        };
    }

}
