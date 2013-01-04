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

package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.datasources.reads.Shard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.iterators.StingSAMIterator;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Collection;

/**
 * @author Joel Thibault
 */
public class ActiveRegionShardDataProvider extends ShardDataProvider {
    final private ReadShardDataProvider readProvider;
    final private LocusShardDataProvider locusProvider;

    public ActiveRegionShardDataProvider(Shard shard, ReadProperties sourceInfo, GenomeLocParser genomeLocParser, StingSAMIterator reads, GenomeLoc locus, LocusIterator locusIterator, IndexedFastaSequenceFile reference, Collection<ReferenceOrderedDataSource> rods) {
        super(shard, genomeLocParser, reference, rods);     // TODO: necessary?
        readProvider = new ReadShardDataProvider(shard, genomeLocParser, reads, reference, rods);
        locusProvider = new LocusShardDataProvider(shard, sourceInfo, genomeLocParser, locus,  locusIterator, reference, rods);
    }

    public ReadShardDataProvider getReadShardDataProvider() {
        return readProvider;
    }

    public LocusShardDataProvider getLocusShardDataProvider(LocusIterator iterator) {
        return locusProvider;
    }
}
