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

import htsjdk.samtools.util.PeekableIterator;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * a ROD view that allows for requests for RODs that overlap intervals on the genome to produce a RefMetaDataTracker
 */
public class IntervalReferenceOrderedView implements ReferenceOrderedView {
    /** a list of the RMDDataState (location->iterators) */
    private final List<RMDDataState> states = new ArrayList<>(1);

    /**
     * Used to get genome locs for reads
     */
    protected final GenomeLocParser genomeLocParser;

    /**
     * The total extent of all reads in this span.  We create iterators from our RODs
     * from the start of this span, to the end.
     */
    private final GenomeLoc shardSpan;

    /**
     * Create a new IntervalReferenceOrderedView taking data from provider and capable of
     * servicing ROD overlap requests within the genomic interval span
     *
     * @param provider a ShardDataProvider to give us data
     * @param span a GenomeLoc span, or null indicating take the entire genome
     */
    public IntervalReferenceOrderedView(final ShardDataProvider provider, final GenomeLoc span) {
        if ( provider == null ) throw new IllegalArgumentException("provider cannot be null");
        if ( provider.hasReferenceOrderedData() && span == null ) throw new IllegalArgumentException("span cannot be null when provider has reference ordered data");

        this.genomeLocParser = provider.getGenomeLocParser();
        this.shardSpan = span;
        provider.register(this);

        // conditional to optimize the case where we don't have any ROD data
        if ( provider.hasReferenceOrderedData() && ! shardSpan.isUnmapped() ) {
            for (final ReferenceOrderedDataSource dataSource : provider.getReferenceOrderedData())
                states.add(new RMDDataState(dataSource, dataSource.seek(shardSpan)));
        }
    }

    /**
     * Testing constructor
     */
    protected IntervalReferenceOrderedView(final GenomeLocParser genomeLocParser,
                                           final GenomeLoc shardSpan,
                                           final List<String> names,
                                           final List<PeekableIterator<RODRecordList>> featureSources) {
        this.genomeLocParser = genomeLocParser;
        this.shardSpan = shardSpan;
        for ( int i = 0; i < names.size(); i++ )
            states.add(new RMDDataState(names.get(i), featureSources.get(i)));
    }

    public Collection<Class<? extends View>> getConflictingViews() {
        List<Class<? extends View>> classes = new ArrayList<>();
        classes.add(ManagingReferenceOrderedView.class);
        return classes;
    }

    /**
     * Get a RefMetaDataTracker containing bindings for all RODs overlapping the start position of loc
     * @param loc a GenomeLoc of size == 1
     * @return a non-null RefMetaDataTracker
     */
    @Override
    public RefMetaDataTracker getReferenceOrderedDataAtLocus(GenomeLoc loc) {
        if ( loc == null ) throw new IllegalArgumentException("loc cannot be null");
        if ( loc.size() != 1 ) throw new IllegalArgumentException("GenomeLoc must have size == 1 but got " + loc);
        return getReferenceOrderedDataForInterval(loc);
    }

    /**
     * Get a RefMetaDataTracker containing bindings for all RODs overlapping interval
     *
     * @param interval a non=null interval
     * @return a non-null RefMetaDataTracker
     */
    public RefMetaDataTracker getReferenceOrderedDataForInterval(final GenomeLoc interval) {
        if ( interval == null ) throw new IllegalArgumentException("Interval cannot be null");

        if ( states.isEmpty() || shardSpan.isUnmapped() ) // optimization for no bindings (common for read walkers)
            return RefMetaDataTracker.EMPTY_TRACKER;
        else {
            final List<RODRecordList> bindings = new ArrayList<>(states.size());
            for ( final RMDDataState state : states )
                bindings.add(state.stream.getOverlapping(interval));
            return new RefMetaDataTracker(bindings);
        }
    }

    /**
     * Trim down all of the ROD managers so that they only hold ROD bindings wit start >= startOfDataToKeep.getStart()
     *
     * @param startOfDataToKeep a non-null genome loc
     */
    public void trimCurrentFeaturesToLoc(final GenomeLoc startOfDataToKeep) {
        if ( startOfDataToKeep == null ) throw new IllegalArgumentException("startOfDataToKeep cannot be null");

        for ( final RMDDataState state : states )
            state.stream.trimCurrentFeaturesToLoc(startOfDataToKeep);
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

    /**
     * Models the traversal state of a given ROD lane.
     */
    private static class RMDDataState {
        public final ReferenceOrderedDataSource dataSource;
        public final IntervalOverlappingRODsFromStream stream;
        private final LocationAwareSeekableRODIterator iterator;

        public RMDDataState(ReferenceOrderedDataSource dataSource, LocationAwareSeekableRODIterator iterator) {
            this.dataSource = dataSource;
            this.iterator = iterator;
            this.stream = new IntervalOverlappingRODsFromStream(dataSource.getName(), new PeekableIterator<>(iterator));
        }

        /**
         * For testing
         */
        public RMDDataState(final String name, final PeekableIterator<RODRecordList> iterator) {
            this.dataSource = null;
            this.iterator = null;
            this.stream = new IntervalOverlappingRODsFromStream(name, new PeekableIterator<>(iterator));
        }

        public void close() {
            if ( dataSource != null )
                dataSource.close( iterator );
        }
    }
}

