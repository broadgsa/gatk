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

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.collections.RODMergingIterator;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileupImpl;

import java.util.*;

/**
 * A view into the reference-ordered data in the provider.
 */
public class RodLocusView extends LocusView implements ReferenceOrderedView {
    /**
     * The data sources along with their current states.
     */
    private RODMergingIterator rodQueue = null;

    Collection<RODRecordList> allTracksHere;

    GenomeLoc lastLoc = null;
    RODRecordList interval = null;

    /**
     * The data sources along with their current states.
     */
    private List<ReferenceOrderedDataState> states = new ArrayList<ReferenceOrderedDataState>();    

    /**
     * Enable debugging output -- todo remove me
     */
    final static boolean DEBUG = false;

    final static String INTERVAL_ROD_NAME = "interval";

    /**
     * Create a new view of reference-ordered data.
     *
     * @param provider
     */
    public RodLocusView( LocusShardDataProvider provider ) {
        super(provider);

        GenomeLoc loc = provider.getLocus();

        List< Iterator<RODRecordList> > iterators = new LinkedList< Iterator<RODRecordList> >();
        for( ReferenceOrderedDataSource dataSource: provider.getReferenceOrderedData() ) {
            if ( DEBUG ) System.out.printf("Shard is %s%n", provider.getLocus());

            // grab the ROD iterator from the data source, and compute the first location in this shard, forwarding
            // the iterator to immediately before it, so that it can be added to the merging iterator primed for
            // next() to return the first real ROD in this shard
            LocationAwareSeekableRODIterator it = dataSource.seek(provider.getLocus());
            it.seekForward(genomeLocParser.createGenomeLoc(loc.getContig(), loc.getStart()-1));

            states.add(new ReferenceOrderedDataState(dataSource,it));            

            // we need to special case the interval so we don't always think there's a rod at the first location
            if ( dataSource.getName().equals(INTERVAL_ROD_NAME) ) {
                if ( interval != null )
                    throw new RuntimeException("BUG: interval local variable already assigned " + interval);
                interval = it.next();
            } else {
                iterators.add( it );
            }
        }

        rodQueue = new RODMergingIterator(iterators);
    }

    @Override
    public RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc ) {
        // special case the interval again -- add it into the ROD
        if ( interval != null ) { allTracksHere.add(interval); }
        return new RefMetaDataTracker(allTracksHere);
    }

    public boolean hasNext() {
        if ( ! rodQueue.hasNext() )
            return false;
        else {
            return ! rodQueue.peekLocation().isPast(locus);
        }
    }

    /**
     * Returns the next covered locus context in the shard.
     * @return Next covered locus context in the shard.
     * @throw NoSuchElementException if no such element exists.
     */
    public AlignmentContext next() {
        if ( DEBUG ) System.out.printf("In RodLocusView.next()...%n");
        RODRecordList datum = rodQueue.next();
        if ( DEBUG ) System.out.printf("In RodLocusView.next(); datum = %s...%n", datum.getLocation());

        if ( DEBUG ) System.out.printf("In RodLocusView.next(): creating tracker...%n");

        allTracksHere = getSpanningTracks(datum);
        GenomeLoc rodSite = datum.getLocation();
        GenomeLoc site = genomeLocParser.createGenomeLoc( rodSite.getContig(), rodSite.getStart(), rodSite.getStart());

        if ( DEBUG ) System.out.printf("rodLocusView.next() is at %s%n", site);

        // calculate the number of skipped bases, and update lastLoc so we can do that again in the next()
        long skippedBases = getSkippedBases( rodSite );
        lastLoc = site;
        return new AlignmentContext(site, new ReadBackedPileupImpl(site), skippedBases);
    }

    private Collection<RODRecordList> getSpanningTracks(RODRecordList marker) {
        return rodQueue.allElementsLTE(marker);
    }

    /**
     * Returns the number of reference bases that have been skipped:
     *
     * 1 -- since the last processed location if we have one
     * 2 -- from the beginning of the shard if this is the first loc
     * 3 -- from the last location to the current position
     *
     * @param currentPos
     * @return
     */
    private long getSkippedBases( GenomeLoc currentPos ) {
        // the minus - is because if lastLoc == null, you haven't yet seen anything in this interval, so it should also be counted as skipped
        Integer compStop = lastLoc == null ? locus.getStart() - 1 : lastLoc.getStop();
        long skippedBases = currentPos.getStart() - compStop  - 1;

        if ( skippedBases < -1 ) { // minus 1 value is ok
            throw new RuntimeException(String.format("BUG: skipped bases=%d is < 0: cur=%s vs. last=%s, shard=%s",
                    skippedBases, currentPos, lastLoc, locus));
        }
        return Math.max(skippedBases, 0);
    }

    /**
     * Get the location one after the last position we will traverse through
     * @return
     */
    public GenomeLoc getLocOneBeyondShard() {
        return genomeLocParser.createGenomeLoc(locus.getContig(),locus.getStop()+1);
    }

    /**
     * How many bases are we skipping from the current location to the end of the interval / shard
     * if we have no more elements
     *
     * @return
     */
    public long getLastSkippedBases() {
        if ( hasNext() )
            throw new RuntimeException("BUG: getLastSkippedBases called when there are elements remaining.");

        return getSkippedBases(getLocOneBeyondShard());
    }

    /**
     * Closes the current view.
     */
    public void close() {
        for( ReferenceOrderedDataState state: states )
            state.dataSource.close( state.iterator );

        rodQueue = null;
        allTracksHere = null;
    }
}