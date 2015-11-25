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

package org.broadinstitute.gatk.engine.datasources.providers;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.engine.ReadProperties;
import org.broadinstitute.gatk.engine.datasources.reads.Shard;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.locusiterator.LocusIterator;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import java.util.Collection;

/**
 * Presents data sharded by locus to the traversal engine.
 *
 * @author mhanna
 * @version 0.1
 */
public class LocusShardDataProvider extends ShardDataProvider {
    /**
     * Information about the source of the read data.
     */
    private final ReadProperties sourceInfo;

    /**
     * The particular locus for which data is provided.  Should be contained within shard.getGenomeLocs().
     */
    private final GenomeLoc locus;

    /**
     * The raw collection of reads.
     */
    private final LocusIterator locusIterator;

    /**
     * Create a data provider for the shard given the reads and reference.
     * @param shard The chunk of data over which traversals happen.
     * @param reference A getter for a section of the reference.
     */
    public LocusShardDataProvider(Shard shard, ReadProperties sourceInfo, GenomeLocParser genomeLocParser, GenomeLoc locus, LocusIterator locusIterator, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        super(shard,genomeLocParser,reference,rods);
        this.sourceInfo = sourceInfo;
        this.locus = locus;
        this.locusIterator = locusIterator;
    }

    /**
     * Returns information about the source of the reads.
     * @return Info about the source of the reads.
     */
    public ReadProperties getSourceInfo() {
        return sourceInfo;
    }

    /**
     * Gets the locus associated with this shard data provider.
     * @return The locus.
     */
    public GenomeLoc getLocus() {
        return locus;
    }

    /**
     * Gets an iterator over all the reads bound by this shard.
     * @return An iterator over all reads in this shard.
     */
    public LocusIterator getLocusIterator() {
        return locusIterator;
    }        

    @Override
    public void close() {
        super.close();
    }
}
