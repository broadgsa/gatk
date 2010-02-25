package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.utils.FlashBackIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MergingIterator;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.iterators.LocusOverflowTracker;

import java.util.*;

/**
 * User: hanna
 * Date: May 21, 2009
 * Time: 2:49:17 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A view into the reference-ordered data in the provider.
 */
public class RodLocusView extends LocusView implements ReferenceOrderedView {
    /**
     * The data sources along with their current states.
     */
    private MergingIterator rodQueue = null;

    RefMetaDataTracker tracker = null;
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
    public RodLocusView( ShardDataProvider provider ) {
        super(provider);

        GenomeLoc loc = provider.getLocus();

        List< Iterator<RODRecordList> > iterators = new LinkedList< Iterator<RODRecordList> >();
        for( ReferenceOrderedDataSource dataSource: provider.getReferenceOrderedData() ) {
            if ( DEBUG ) System.out.printf("Shard is %s%n", provider.getLocus());

            // grab the ROD iterator from the data source, and compute the first location in this shard, forwarding
            // the iterator to immediately before it, so that it can be added to the merging iterator primed for
            // next() to return the first real ROD in this shard
            FlashBackIterator it = (FlashBackIterator)dataSource.seek(provider.getShard());
            it.seekForward(GenomeLocParser.createGenomeLoc(loc.getContigIndex(), loc.getStart()-1));

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

        rodQueue = new MergingIterator(iterators);

        //throw new StingException("RodLocusView currently disabled");
    }

    public RefMetaDataTracker getReferenceOrderedDataAtLocus( GenomeLoc loc ) {
        return tracker;
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

        // Update the tracker here for use
        Collection<RODRecordList> allTracksHere = getSpanningTracks(datum);
        tracker = createTracker(allTracksHere);

        GenomeLoc rodSite = datum.getLocation();
        GenomeLoc site = GenomeLocParser.createGenomeLoc( rodSite.getContigIndex(), rodSite.getStart(), rodSite.getStart());

        if ( DEBUG ) System.out.printf("rodLocusView.next() is at %s%n", site);

        // calculate the number of skipped bases, and update lastLoc so we can do that again in the next()
        long skippedBases = getSkippedBases( rodSite );
        lastLoc = site;
        return new AlignmentContext(site, new ReadBackedPileup(site), skippedBases);
    }

    public LocusOverflowTracker getLocusOverflowTracker() {
        // we don't have an overflow tracker
        return null;  
    }

    private RefMetaDataTracker createTracker( Collection<RODRecordList> allTracksHere ) {
        RefMetaDataTracker t = new RefMetaDataTracker();
        for ( RODRecordList track : allTracksHere ) {
            if ( ! t.hasROD(track.getName()) )
                t.bind(track.getName(), track);
        }

        // special case the interval again -- add it into the ROD
        if ( interval != null ) { t.bind(interval.getName(), interval); }

        return t;
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
        Long compStop = lastLoc == null ? locus.getStart() - 1 : lastLoc.getStop();
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
        return GenomeLocParser.createGenomeLoc(locus.getContigIndex(),locus.getStop()+1);
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

    public RefMetaDataTracker getTracker() {
        return tracker;
    }

    /**
     * Closes the current view.
     */
    public void close() {
        for( ReferenceOrderedDataState state: states )
            state.dataSource.close( state.iterator );

        rodQueue = null;
        tracker = null;
    }
}