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

package org.broadinstitute.gatk.engine.walkers;

import htsjdk.samtools.SAMRecord;

import java.util.Collection;

/**
 * Walks over all pairs/collections of reads in a BAM file sorted by
 * read name.
 *
 * @author mhanna
 * @version 0.1
 */
@Requires({DataSource.READS})
public abstract class ReadPairWalker<MapType,ReduceType> extends Walker<MapType,ReduceType> {

    /**
     * Optionally filters out read pairs.
     * @param reads collections of all reads with the same read name.
     * @return True to process the reads with map/reduce; false otherwise.
     */
    public boolean filter(Collection<SAMRecord> reads) {
        // Keep all pairs by default.
        return true;
    }

    /**
     * Maps a read pair to a given reduce of type MapType.  Semantics determined by subclasser.
     * @param reads Collection of reads having the same name.
     * @return Semantics defined by implementer.
     */
    public abstract MapType map(Collection<SAMRecord> reads);

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);

}
