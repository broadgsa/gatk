/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.providers;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.datasources.reads.ReadShard;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/** a ROD view for reads. This provides the Read traversals a way of getting a RefMetaDataTracker */
public class ReadBasedReferenceOrderedView implements View {
    // a list of the RMDDataState (location->iterators)
    private final List<RMDDataState> states = new ArrayList<RMDDataState>(1);
    private final static RefMetaDataTracker EMPTY_TRACKER = new RefMetaDataTracker();

    /**
     * Used to get genome locs for reads
     */
    private final GenomeLocParser genomeLocParser;

    /**
     * The total extent of all reads in this span.  We create iterators from our RODs
     * from the start of this span, to the end.
     */
    private final GenomeLoc shardSpan;

    public ReadBasedReferenceOrderedView(final ShardDataProvider provider) {
        this.genomeLocParser = provider.getGenomeLocParser();
        // conditional to optimize the case where we don't have any ROD data
        this.shardSpan = provider.getReferenceOrderedData() != null ? ((ReadShard)provider.getShard()).getReadsSpan() : null;
        provider.register(this);

        if ( provider.getReferenceOrderedData() != null && ! shardSpan.isUnmapped() ) {
            for (ReferenceOrderedDataSource dataSource : provider.getReferenceOrderedData())
                states.add(new RMDDataState(dataSource, dataSource.seek(shardSpan)));
        }
    }


    /**
     * Testing constructor
     */
    protected ReadBasedReferenceOrderedView(final GenomeLocParser genomeLocParser,
                                            final GenomeLoc shardSpan,
                                            final List<String> names,
                                            final List<PeekableIterator<RODRecordList>> featureSources) {
        this.genomeLocParser = genomeLocParser;
        this.shardSpan = shardSpan;
        for ( int i = 0; i < names.size(); i++ )
            states.add(new RMDDataState(names.get(i), featureSources.get(i)));
    }

    public Collection<Class<? extends View>> getConflictingViews() {
        List<Class<? extends View>> classes = new ArrayList<Class<? extends View>>();
        classes.add(ManagingReferenceOrderedView.class);
        return classes;
    }

    /**
     * create a RefMetaDataTracker given the current read
     *
     * @param rec the read
     *
     * @return a RefMetaDataTracker for the read, from which you can get ROD -> read alignments
     */
    @Requires("rec != null")
    @Ensures("result != null")
    public RefMetaDataTracker getReferenceOrderedDataForRead(final SAMRecord rec) {
        if ( rec.getReadUnmappedFlag() )
            // empty RODs for unmapped reads
            return new RefMetaDataTracker();
        else
            return getReferenceOrderedDataForInterval(genomeLocParser.createGenomeLoc(rec));
    }

    @Requires({"interval != null", "shardSpan == null || shardSpan.isUnmapped() || shardSpan.containsP(interval)"})
    @Ensures("result != null")
    public RefMetaDataTracker getReferenceOrderedDataForInterval(final GenomeLoc interval) {
        if ( states.isEmpty() || shardSpan.isUnmapped() ) // optimization for no bindings (common for read walkers)
            return EMPTY_TRACKER;
        else {
            final List<RODRecordList> bindings = new ArrayList<RODRecordList>(states.size());
            for ( final RMDDataState state : states )
                bindings.add(state.stream.getOverlapping(interval));
            return new RefMetaDataTracker(bindings);
        }
    }

    /**
     * Closes the current view.
     */
    public void close() {
        for (final RMDDataState state : states)
            state.close();

        // Clear out the existing data so that post-close() accesses to this data will fail-fast.
        states.clear();
    }

    /** Models the traversal state of a given ROD lane. */
    private static class RMDDataState {
        public final ReferenceOrderedDataSource dataSource;
        public final IntervalOverlappingRODsFromStream stream;
        private final LocationAwareSeekableRODIterator iterator;

        public RMDDataState(ReferenceOrderedDataSource dataSource, LocationAwareSeekableRODIterator iterator) {
            this.dataSource = dataSource;
            this.iterator = iterator;
            this.stream = new IntervalOverlappingRODsFromStream(dataSource.getName(), new PeekableIterator<RODRecordList>(iterator));
        }

        /**
         * For testing
         */
        public RMDDataState(final String name, final PeekableIterator<RODRecordList> iterator) {
            this.dataSource = null;
            this.iterator = null;
            this.stream = new IntervalOverlappingRODsFromStream(name, new PeekableIterator<RODRecordList>(iterator));
        }

        public void close() {
            if ( dataSource != null )
                dataSource.close( iterator );
        }
    }
}

