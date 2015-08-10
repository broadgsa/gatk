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

/**
 * Batch granular file pointers into potentially larger shards.
 */
public class LocusShardBalancer extends ShardBalancer {
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
                FilePointer current = filePointers.next();

                // FilePointers have already been combined as necessary at the IntervalSharder level. No
                // need to do so again here.

                return new LocusShard(parser,readsDataSource,current.getLocations(),current.fileSpans);
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from shard balancing iterator");
            }
        };
    }
}
