/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.traversals;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionTraversalParameters;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.activeregion.*;
import org.broadinstitute.sting.utils.progressmeter.ProgressMeter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.*;

/**
 * Implement active region traversal
 *
 * User: depristo
 * Date: 1/9/13
 * Time: 4:45 PM
 *
 * Live region:
 *
 *   The ART tracks a thing called the live region.  The live region is a position on a specific contig
 *   of the alignment start of the last read we processed during this traversal.  Because the
 *   read stream is sorted, future reads must occurs in the the live region.  Therefore the the dead region
 *   (everything to the left of the live boundary) cannot have any more read data.  The live / dead
 *   regions are used to decide when we can safely call map on active regions, as only active regions
 *   contained completely within the dead region (including extensions) have a complete set of read data
 *   in the collected read list.  All of the data related to the live region is captured by the local
 *   variable spanOfLastReadSeen
 *
 */
public class TraverseActiveRegions<M, T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,LocusShardDataProvider> {
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);
    protected final static boolean LOG_READ_CARRYING = false;

    // set by the traversal
    private boolean walkerHasPresetRegions = false;
    private int activeRegionExtension = -1;
    private int maxRegionSize = -1;
    private int minRegionSize = -1;

    private final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();

    private TAROrderedReadCache myReads = null;

    private GenomeLoc spanOfLastReadSeen = null;
    private ActivityProfile activityProfile = null;
    int maxReadsInMemory = 0;
    ActiveRegionWalker walker;

    /**
     * Have the debugging output streams been initialized already?
     *
     * We have to do lazy initialization because when the initialize() function is called
     * the streams aren't yet initialized in the GATK walker.
     */
    private boolean streamsInitialized = false;

    @Override
    public void initialize(GenomeAnalysisEngine engine, Walker walker, ProgressMeter progressMeter) {
        super.initialize(engine, walker, progressMeter);

        this.walker = (ActiveRegionWalker)walker;
        if ( this.walker.wantsExtendedReads() && ! this.walker.wantsNonPrimaryReads() ) {
            throw new IllegalArgumentException("Active region walker " + this.walker + " requested extended events but not " +
                    "non-primary reads, an inconsistent state.  Please modify the walker");
        }

        ActiveRegionTraversalParameters annotation = walker.getClass().getAnnotation(ActiveRegionTraversalParameters.class);
        this.activeRegionExtension = this.walker.activeRegionExtension == null ? annotation.extension() : this.walker.activeRegionExtension;
        this.maxRegionSize = this.walker.activeRegionMaxSize == null ? annotation.maxRegion() : this.walker.activeRegionMaxSize;
        this.minRegionSize = annotation.minRegion();
        final double bandPassSigma = this.walker.bandPassSigma == null ? annotation.bandPassSigma() : this.walker.bandPassSigma;
        walkerHasPresetRegions = this.walker.hasPresetActiveRegions();

        activityProfile = new BandPassActivityProfile(engine.getGenomeLocParser(), engine.getIntervals(), BandPassActivityProfile.MAX_FILTER_SIZE, bandPassSigma);

        if ( walkerHasPresetRegions ) {
            // we load all of the preset locations into the workQueue
            for ( final GenomeLoc loc : this.walker.getPresetActiveRegions()) {
                workQueue.add(new ActiveRegion(loc, null, true, engine.getGenomeLocParser(), getActiveRegionExtension()));
            }
        }

        final int maxReadsAcrossSamples = annotation.maxReadsToHoldInMemoryPerSample() * SampleUtils.getSAMFileSamples(engine).size();
        final int maxReadsToHoldInMemory = Math.min(maxReadsAcrossSamples, annotation.maxReadsToHoldTotal());
        myReads = new TAROrderedReadCache(maxReadsToHoldInMemory);
    }

    // -------------------------------------------------------------------------------------
    //
    // Utility functions
    //
    // -------------------------------------------------------------------------------------

    protected int getActiveRegionExtension() {
        return activeRegionExtension;
    }

    protected int getMaxRegionSize() {
        return maxRegionSize;
    }

    protected int getMinRegionSize() {
        return minRegionSize;
    }

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public String toString() {
        return "TraverseActiveRegions";
    }

    /**
     * Is the loc outside of the intervals being requested for processing by the GATK?
     * @param loc
     * @return
     */
    protected boolean outsideEngineIntervals(final GenomeLoc loc) {
        return engine.getIntervals() != null && ! engine.getIntervals().overlaps(loc);
    }

    protected ReferenceOrderedView getReferenceOrderedView(final ActiveRegionWalker<M, T> walker,
                                                           final LocusShardDataProvider dataProvider,
                                                           final LocusView locusView) {
        if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
            return new ManagingReferenceOrderedView( dataProvider );
        else
            return (RodLocusView)locusView;
    }


    // -------------------------------------------------------------------------------------
    //
    // Actual traverse function
    //
    // -------------------------------------------------------------------------------------

    /**
     * Did read appear in the last shard?
     *
     * When we transition across shard boundaries we see duplicate reads because
     * each shard contains the reads that *overlap* the shard.  So if we just finished
     * shard 1-1000 and are now in 1001-2000 we'll see duplicate reads from 1001
     * that overlapped 1-1000.  This function tests read to determine if we would have
     * seen it before by asking if read.getAlignmentStart() is less than the
     * stop position of the last seen read at the start of the traversal.  The reason
     * we need to use the location of the last read at the start of the traversal
     * is that we update the lastRead during the traversal, and we only want to filter
     * out reads whose start is before the last read of the previous shard, not the
     * current shard.
     *
     * @param locOfLastReadAtTraversalStart the location of the last read seen at the start of the traversal
     * @param read the read we want to test if it's already been seen in the last shard
     * @return true if read would have appeared in the last shard, false otherwise
     */
    @Requires({"read != null"})
    private boolean appearedInLastShard(final GenomeLoc locOfLastReadAtTraversalStart, final GATKSAMRecord read) {
        if ( locOfLastReadAtTraversalStart == null )
            // we're in the first shard, so obviously the answer is no
            return false;
        else {
            // otherwise check to see if the alignment occurred in the previous shard
            return read.getAlignmentStart() <= locOfLastReadAtTraversalStart.getStart()
                    // we're on the same contig
                    && read.getReferenceIndex() == locOfLastReadAtTraversalStart.getContigIndex();
        }

    }

    @Override
    public T traverse( final ActiveRegionWalker<M,T> walker,
                       final LocusShardDataProvider dataProvider,
                       T sum) {
        if ( LOG_READ_CARRYING || logger.isDebugEnabled() )
            logger.info(String.format("TraverseActiveRegions.traverse: Shard is %s", dataProvider));

        final LocusView locusView = new AllLocusView(dataProvider);
        final LocusReferenceView referenceView = new LocusReferenceView( walker, dataProvider );
        final ReferenceOrderedView referenceOrderedDataView = getReferenceOrderedView(walker, dataProvider, locusView);

        // We keep processing while the next reference location is within the interval
        final GenomeLoc locOfLastReadAtTraversalStart = spanOfLastSeenRead();

        while( locusView.hasNext() ) {
            final AlignmentContext locus = locusView.next();
            final GenomeLoc location = locus.getLocation();

            rememberLastLocusLocation(location);

            // get all of the new reads that appear in the current pileup, and them to our list of reads
            // provided we haven't seen them before
            final Collection<GATKSAMRecord> reads = locusView.getLIBS().transferReadsFromAllPreviousPileups();
            for( final GATKSAMRecord read : reads ) {
                if ( ! appearedInLastShard(locOfLastReadAtTraversalStart, read) ) {
                    rememberLastReadLocation(read);
                    myReads.add(read);
                }
            }

            // skip this location -- it's not part of our engine intervals
            if ( outsideEngineIntervals(location) )
                continue;

            // we've move across some interval boundary, restart profile
            final boolean flushProfile = ! activityProfile.isEmpty()
                    && ( activityProfile.getContigIndex() != location.getContigIndex()
                       || location.getStart() != activityProfile.getStop() + 1);
            sum = processActiveRegions(walker, sum, flushProfile, false);

            dataProvider.getShard().getReadMetrics().incrementNumIterations();

            // create reference context. Note that if we have a pileup of "extended events", the context will
            // hold the (longest) stretch of deleted reference bases (if deletions are present in the pileup).
            final ReferenceContext refContext = referenceView.getReferenceContext(location);

            // Iterate forward to get all reference ordered data covering this location
            final RefMetaDataTracker tracker = referenceOrderedDataView.getReferenceOrderedDataAtLocus(locus.getLocation(), refContext);

            // Call the walkers isActive function for this locus and add them to the list to be integrated later
            addIsActiveResult(walker, tracker, refContext, locus);

            maxReadsInMemory = Math.max(myReads.size(), maxReadsInMemory);
            printProgress(location);
        }

        updateCumulativeMetrics(dataProvider.getShard());

        return sum;
    }

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public T endTraversal(final Walker<M, T> walker, T sum) {
        return processActiveRegions((ActiveRegionWalker<M, T>)walker, sum, true, true);
    }

    // -------------------------------------------------------------------------------------
    //
    // Functions to manage and interact with the live / dead zone
    //
    // -------------------------------------------------------------------------------------

    /**
     * Update the live region to reflect that the last read we've seen in the traversal is read
     *
     * Requires that sequential calls always be provided reads in coordinate sorted order
     *
     * @param read the last read we've seen during the traversal
     */
    @Requires({"read != null"})
    protected void rememberLastReadLocation(final GATKSAMRecord read) {
        final GenomeLoc currentLocation = engine.getGenomeLocParser().createGenomeLoc(read);
        if ( spanOfLastReadSeen == null )
            spanOfLastReadSeen = currentLocation;
        else {
            if ( currentLocation.isBefore(spanOfLastReadSeen) )
                throw new IllegalStateException("Updating last read seen in the traversal with read " + read + " with span " + currentLocation + " but this occurs before the previously seen read " + spanOfLastReadSeen);
            spanOfLastReadSeen = currentLocation;
        }
    }

    /**
     * Update the live region to reflect that we've reached locus
     *
     * This function is complementary to #rememberLastReadLocation, but if we don't have any reads for a long
     * time (e.g., there's no coverage) we will keep active regions around far longer than necessary.
     *
     * Only updates the span if it's beyond the last seen
     *
     * @param currentLocation the current location we've processed on the genome
     */
    protected void rememberLastLocusLocation(final GenomeLoc currentLocation) {
        if ( spanOfLastReadSeen == null )
            spanOfLastReadSeen = currentLocation;
        else {
            if ( currentLocation.isPast(spanOfLastReadSeen) )
                spanOfLastReadSeen = currentLocation;
        }
    }


    /**
     * Get a GenomeLoc indicating the start (heading to the right) of the live ART region.
     * @return the left-most position of the live region on the genome
     */
    protected GenomeLoc spanOfLastSeenRead() {
        return spanOfLastReadSeen;
    }

    /**
     * Is the active region completely within the traversal's dead zone?
     *
     * @param region the region we want to test
     * @return true if the extended location of region is completely within the current dead zone, false otherwise
     */
    protected boolean regionCompletelyWithinDeadZone(final ActiveRegion region) {
        if ( spanOfLastSeenRead() == null )
            return false;

        final int contigCmp = region.getExtendedLoc().compareContigs(spanOfLastSeenRead());
        if ( contigCmp > 0 )
            throw new IllegalStateException("Active region " + region + " on a contig after last seen read " + spanOfLastSeenRead());
        else {
            return contigCmp < 0 || region.getExtendedLoc().getStop() < spanOfLastSeenRead().getStart();
        }
    }

    /**
     * Is the read dead?  That is, can it no longer be in any future active region, and therefore can be discarded?
     *
     * read: start |--------> stop ------ stop + extension
     * region:                      start |-----------------| end
     *
     * Since the regions are coming in order, read could potentially be contained in a future interval if
     * stop + activeRegionExtension >= end.  If, on the other hand, stop + extension is < the end
     * of this region, then we can discard it, since any future region could only include reads
     * up to end + 1 - extension.
     *
     * Note that this function doesn't care about the dead zone.  We're assuming that by
     * actually calling this function with an active region that region is already in the dead zone,
     * so checking that the read is in the dead zone doesn't make sense.
     *
     * @param read the read we're testing
     * @param activeRegion the current active region
     * @return true if the read is dead, false other
     */
    @Requires({"read != null", "activeRegion != null"})
    private boolean readCannotOccurInAnyMoreActiveRegions(final GATKSAMRecord read, final ActiveRegion activeRegion) {
        return read.getReferenceIndex() < activeRegion.getLocation().getContigIndex() ||
                ( read.getReferenceIndex() == activeRegion.getLocation().getContigIndex()
                        && read.getAlignmentEnd() + getActiveRegionExtension() < activeRegion.getLocation().getStop() );
    }

    // -------------------------------------------------------------------------------------
    //
    // Functions to write out activity profiles and active regions
    //
    // -------------------------------------------------------------------------------------

    /**
     * Initialize the debugging output streams (activity profile and active regions), if not done so already
     */
    @Ensures("streamsInitialized == true")
    private void initializeOutputStreamsIfNecessary() {
        if ( ! streamsInitialized ) {
            streamsInitialized = true;
            if ( walker.activityProfileOutStream != null ) {
                printIGVFormatHeader(walker.activityProfileOutStream, "line", "ActivityProfile");
            }

            if ( walker.activeRegionOutStream != null ) {
                printIGVFormatHeader(walker.activeRegionOutStream, "line", "ActiveRegions");
            }
        }
    }

    /**
     * Helper function to write out a IGV formatted line to out, at loc, with values
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * @param out a non-null PrintStream where we'll write our line
     * @param graphType the type of graph to show in IGV for this track
     * @param columns the column names for this IGV track
     */
    @Requires({
            "out != null",
            "graphType != null",
            "columns.length > 0"
    })
    private void printIGVFormatHeader(final PrintStream out, final String graphType, final String ... columns ) {
        out.printf("#track graphType=%s%n", graphType);
        out.printf("Chromosome\tStart\tEnd\tFeature\t%s%n", Utils.join("\t", columns));

    }

    /**
     * Helper function to write out a IGV formatted line to out, at loc, with values
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * @param out a non-null PrintStream where we'll write our line
     * @param loc the location of values
     * @param featureName string name of this feature (see IGV format)
     * @param values the floating point values to associate with loc and feature name in out
     */
    @Requires({
            "out != null",
            "loc != null",
            "values.length > 0"
    })
    private void printIGVFormatRow(final PrintStream out, final GenomeLoc loc, final String featureName, final double ... values) {
        // note that start and stop are 0 based, but the stop is exclusive so we don't subtract 1
        out.printf("%s\t%d\t%d\t%s", loc.getContig(), loc.getStart() - 1, loc.getStop(), featureName);
        for ( final double value : values )
            out.print(String.format("\t%.3f", value));
        out.println();
    }

    /**
     * Write out activity profile information, if requested by the walker
     *
     * @param states the states in the current activity profile
     */
    @Requires("states != null")
    private void writeActivityProfile(final List<ActivityProfileState> states) {
        if ( walker.activityProfileOutStream != null ) {
            initializeOutputStreamsIfNecessary();
            for ( final ActivityProfileState state : states ) {
                printIGVFormatRow(walker.activityProfileOutStream, state.getLoc(), "state", Math.min(state.isActiveProb, 1.0));
            }
        }
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param region the region we're currently operating on
     */
    @Requires("region != null")
    private void writeActiveRegion(final ActiveRegion region) {
        if( walker.activeRegionOutStream != null ) {
            initializeOutputStreamsIfNecessary();
            printIGVFormatRow(walker.activeRegionOutStream, region.getLocation().getStartLocation(),
                    "end-marker", 0.0);
            printIGVFormatRow(walker.activeRegionOutStream, region.getLocation(),
                    "size=" + region.getLocation().size(), region.isActive() ? 1.0 : -1.0);
        }
    }


    // -------------------------------------------------------------------------------------
    //
    // Functions to process active regions that are ready for map / reduce calls
    //
    // -------------------------------------------------------------------------------------

    /**
     * Invoke the walker isActive function, and incorporate its result into the activity profile
     *
     * @param walker the walker we're running
     * @param tracker the ref meta data tracker to pass on to the isActive function of walker
     * @param refContext the refContext to pass on to the isActive function of walker
     * @param locus the AlignmentContext to pass on to the isActive function of walker
     */
    private void addIsActiveResult(final ActiveRegionWalker<M, T> walker,
                                   final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                   final AlignmentContext locus) {
        // must be called, even if we won't use the result, to satisfy walker contract
        final ActivityProfileState state = walker.isActive( tracker, refContext, locus );
        if ( walker.forceActive) state.isActiveProb = 1.0;
        if ( ! walkerHasPresetRegions ) {
            activityProfile.add(state);
        }
    }

    /**
     * Take the individual isActive calls and integrate them into contiguous active regions and
     * add these blocks of work to the work queue
     * band-pass filter the list of isActive probabilities and turn into active regions
     */
    private T processActiveRegions(final ActiveRegionWalker<M, T> walker, T sum, final boolean flushActivityProfile, final boolean forceAllRegionsToBeActive) {
        if ( ! walkerHasPresetRegions ) {
            // We don't have preset regions, so we get our regions from the activity profile
            final Collection<ActiveRegion> activeRegions = activityProfile.popReadyActiveRegions(getActiveRegionExtension(), getMinRegionSize(), getMaxRegionSize(), flushActivityProfile);
            workQueue.addAll(activeRegions);
            if ( ! activeRegions.isEmpty() && logger.isDebugEnabled() ) logger.debug("Integrated " + activityProfile.size() + " isActive calls into " + activeRegions.size() + " regions." );
        }

        // Since we've traversed sufficiently past this point (or this contig!) in the workQueue we can unload those regions and process them
        while( workQueue.peek() != null ) {
            final ActiveRegion activeRegion = workQueue.peek();
            if ( forceAllRegionsToBeActive || regionCompletelyWithinDeadZone(activeRegion) ) {
                writeActivityProfile(activeRegion.getSupportingStates());
                writeActiveRegion(activeRegion);
                sum = processActiveRegion( workQueue.remove(), sum, walker );
            } else {
                break;
            }
        }

        return sum;
    }

    private T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker) {
        final List<GATKSAMRecord> stillLive = new LinkedList<GATKSAMRecord>();
        for ( final GATKSAMRecord read : myReads.popCurrentReads() ) {
            boolean killed = false;
            final GenomeLoc readLoc = this.engine.getGenomeLocParser().createGenomeLoc( read );

            if( activeRegion.getLocation().overlapsP( readLoc ) ) {
                activeRegion.add(read);

                if ( ! walker.wantsNonPrimaryReads() ) {
                    killed = true;
                }
            } else if( walker.wantsExtendedReads() && activeRegion.getExtendedLoc().overlapsP( readLoc )) {
                activeRegion.add( read );
            }

            // if the read hasn't already been killed, check if it cannot occur in any more active regions, and maybe kill it
            if ( ! killed && readCannotOccurInAnyMoreActiveRegions(read, activeRegion) ) {
                killed = true;
            }

            // keep track of all of the still live active regions
            if ( ! killed ) stillLive.add(read);
        }
        myReads.addAll(stillLive);

        if ( logger.isDebugEnabled() ) {
            logger.debug(">> Map call with " + activeRegion.getReads().size() + " " + (activeRegion.isActive() ? "active" : "inactive") + " reads @ " + activeRegion.getLocation() + " with full extent: " + activeRegion.getReadSpanLoc());
        }

        if ( LOG_READ_CARRYING )
            logger.info(String.format("Processing region %20s span=%3d active?=%5b with %4d reads.  Overall max reads carried is %s",
                    activeRegion.getLocation(), activeRegion.getLocation().size(), activeRegion.isActive(), activeRegion.size(), maxReadsInMemory));

        final M x = walker.map(activeRegion, null);
        return walker.reduce( x, sum );
    }
}
