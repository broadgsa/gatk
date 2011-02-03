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

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeMap;

/** a ROD view for reads. This provides the Read traversals a way of getting a ReadMetaDataTracker */
public class ReadBasedReferenceOrderedView implements View {
    private final WindowedData window;

    public ReadBasedReferenceOrderedView(ShardDataProvider provider) {
        window = new WindowedData(provider);
        provider.register(this);
    }

    /**
     * for testing only please
     *
     * @param data the window provider
     */
    ReadBasedReferenceOrderedView(WindowedData data) {
        window = data;
    }

    public ReadMetaDataTracker getReferenceOrderedDataForRead(SAMRecord read) {
        return window.getTracker(read);
    }

    public Collection<Class<? extends View>> getConflictingViews() {
        List<Class<? extends View>> classes = new ArrayList<Class<? extends View>>();
        classes.add(ManagingReferenceOrderedView.class);
        return classes;
    }

    public void close() {
        if (window != null) window.close();
    }
}


/** stores a window of data, dropping RODs if we've passed the new reads start point. */
class WindowedData {
    // the queue of possibly in-frame RODs; RODs are removed as soon as they are out of scope
    private final TreeMap<Integer, RODMetaDataContainer> mapping = new TreeMap<Integer, RODMetaDataContainer>();

    // our current location from the last read we processed
    private GenomeLoc currentLoc;

    // a list of the RMDDataState (location->iterators)
    private List<RMDDataState> states;

    // the provider; where we get all our information
    private final ShardDataProvider provider;

    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(WindowedData.class);

    /**
     * create a WindowedData given a shard provider
     *
     * @param provider the ShardDataProvider
     */
    public WindowedData(ShardDataProvider provider) {
        this.provider = provider;
    }

    /**
     * load the states dynamically, since the only way to get a genome loc is from the read (the shard doesn't have one)
     *
     * @param provider the ShardDataProvider
     * @param rec      the current read
     */
    private void getStates(ShardDataProvider provider, SAMRecord rec) {

        int stop = Integer.MAX_VALUE;
        // figure out the appropriate alignment stop
        if (provider.hasReference()) {
            stop = provider.getReference().getSequenceDictionary().getSequence(rec.getReferenceIndex()).getSequenceLength();
        }

        // calculate the range of positions we need to look at
        GenomeLoc range = provider.getGenomeLocParser().createGenomeLoc(rec.getReferenceName(),
                rec.getAlignmentStart(),
                stop);
        states = new ArrayList<RMDDataState>();
        if (provider != null && provider.getReferenceOrderedData() != null)
            for (ReferenceOrderedDataSource dataSource : provider.getReferenceOrderedData())
                states.add(new RMDDataState(dataSource, dataSource.seek(range)));
    }

    /**
     * this function is for testing only
     *
     * @param states a  list of RMDDataState to initialize with
     */
    WindowedData(List<RMDDataState> states) {
        this.states = states;
        provider = null;
    }

    /**
     * create a ReadMetaDataTracker given the current read
     *
     * @param rec the read
     *
     * @return a ReadMetaDataTracker for the read, from which you can get ROD -> read alignments
     */
    public ReadMetaDataTracker getTracker(SAMRecord rec) {
        updatePosition(rec);
        return new ReadMetaDataTracker(provider.getGenomeLocParser(), rec, mapping);
    }

    /**
     * update the position we're storing
     *
     * @param rec the read to use for start and end
     */
    private void updatePosition(SAMRecord rec) {
        if (states == null) getStates(this.provider, rec);
        currentLoc = provider.getGenomeLocParser().createGenomeLoc(rec);

        // flush the queue looking for records we've passed over
        while (mapping.size() > 0 && mapping.firstKey() < currentLoc.getStart())
            mapping.pollFirstEntry(); // toss away records that we've passed

        // add new data to the queue
        for (RMDDataState state : states) {
            // move into position
            while (state.iterator.hasNext() && state.iterator.peekNextLocation().isBefore(currentLoc))
                state.iterator.next();
            while (state.iterator.hasNext() && state.iterator.peekNextLocation().overlapsP(currentLoc)) {
                RODRecordList list = state.iterator.next();
                for (GATKFeature datum : list) {
                    if (!mapping.containsKey(list.getLocation().getStart()))
                        mapping.put(list.getLocation().getStart(), new RODMetaDataContainer());
                    mapping.get(list.getLocation().getStart()).addEntry(datum);
                }
            }
        }
    }

    /** Closes the current view. */
    public void close() {
        if (states == null) return;
        for (RMDDataState state : states)
            state.dataSource.close( state.iterator );

        // Clear out the existing data so that post-close() accesses to this data will fail-fast.
        states = null;
    }


}

/** Models the traversal state of a given ROD lane. */
class RMDDataState {
    public final ReferenceOrderedDataSource dataSource;
    public final LocationAwareSeekableRODIterator iterator;

    public RMDDataState(ReferenceOrderedDataSource dataSource, LocationAwareSeekableRODIterator iterator) {
        this.dataSource = dataSource;
        this.iterator = iterator;
    }
}
