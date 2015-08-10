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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

/**
 * ActiveRegionShardBalancer
 *
 * Merges all of the file pointer information for a single contig index into a single
 * combined shard.  The purpose of doing this is to ensure that the HaplotypeCaller, which
 * doesn't support TreeReduction by construction, gets all of the data on a single
 * contig together so the the NanoSchedule runs efficiently
 */
public class ActiveRegionShardBalancer extends ShardBalancer {
    /**
     * Convert iterators of file pointers into balanced iterators of shards.
     * @return An iterator over balanced shards.
     */
    public Iterator<Shard> iterator() {
        return new Iterator<Shard>() {
            public boolean hasNext() {
                return filePointers.hasNext();
            }

            public Shard next() {
                FilePointer current = getCombinedFilePointersOnSingleContig();

                // FilePointers have already been combined as necessary at the IntervalSharder level. No
                // need to do so again here.

                return new LocusShard(parser,readsDataSource,current.getLocations(),current.fileSpans);
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from shard balancing iterator");
            }
        };
    }

    /**
     * Combine all of the file pointers in the filePointers iterator into a single combined
     * FilePointer that spans all of the file pointers on a single contig
     * @return a non-null FilePointer
     */
    private FilePointer getCombinedFilePointersOnSingleContig() {
        FilePointer current = filePointers.next();

        final List<FilePointer> toCombine = new LinkedList<>();
        toCombine.add(current);

        while ( filePointers.hasNext() &&
                current.isRegionUnmapped == filePointers.peek().isRegionUnmapped &&
                (current.getContigIndex() == filePointers.peek().getContigIndex() || current.isRegionUnmapped) ) {
            toCombine.add(filePointers.next());
        }

        return FilePointer.union(toCombine, parser);
    }
}
