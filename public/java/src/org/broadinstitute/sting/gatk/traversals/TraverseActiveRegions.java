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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionExtension;
import org.broadinstitute.sting.gatk.walkers.ActiveRegionWalker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActivityProfile;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.progressmeter.ProgressMeter;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: depristo
 * Date: 1/9/13
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TraverseActiveRegions<M, T> extends TraversalEngine<M,T,ActiveRegionWalker<M,T>,LocusShardDataProvider> {
    protected final static boolean DEBUG = false;

    // set by the tranversal
    private int activeRegionExtension = -1;
    private int maxRegionSize = -1;

    /**
     * our log, which we want to capture anything from this class
     */
    protected final static Logger logger = Logger.getLogger(TraversalEngine.class);
    protected final LinkedList<ActiveRegion> workQueue = new LinkedList<ActiveRegion>();

    abstract protected T processActiveRegion(final ActiveRegion activeRegion, final T sum, final ActiveRegionWalker<M, T> walker);

    /**
     * Special function called in LinearMicroScheduler to empty out the work queue.
     * Ugly for now but will be cleaned up when we push this functionality more into the engine
     */
    public abstract T endTraversal(final Walker<M, T> walker, T sum);

    protected int getActiveRegionExtension() {
        return activeRegionExtension;
    }

    protected int getMaxRegionSize() {
        return maxRegionSize;
    }

    @Override
    public String getTraversalUnits() {
        return "active regions";
    }

    @Override
    public void initialize(GenomeAnalysisEngine engine, Walker walker, ProgressMeter progressMeter) {
        super.initialize(engine, walker, progressMeter);
        activeRegionExtension = walker.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        maxRegionSize = walker.getClass().getAnnotation(ActiveRegionExtension.class).maxRegion();

        final ActiveRegionWalker arWalker = (ActiveRegionWalker)walker;
        if ( arWalker.wantsExtendedReads() && ! arWalker.wantsNonPrimaryReads() ) {
            throw new IllegalArgumentException("Active region walker " + arWalker + " requested extended events but not " +
                    "non-primary reads, an inconsistent state.  Please modify the walker");
        }
    }

    /**
     * Is the loc outside of the intervals being requested for processing by the GATK?
     * @param loc
     * @return
     */
    protected boolean outsideEngineIntervals(final GenomeLoc loc) {
        return engine.getIntervals() != null && ! engine.getIntervals().overlaps(loc);
    }

    /**
     * Take the individual isActive calls and integrate them into contiguous active regions and
     * add these blocks of work to the work queue
     * band-pass filter the list of isActive probabilities and turn into active regions
     *
     * @param profile
     * @param activeRegions
     * @return
     */
    protected ActivityProfile incorporateActiveRegions(final ActivityProfile profile,
                                                       final List<ActiveRegion> activeRegions) {
        if ( profile.isEmpty() )
            throw new IllegalStateException("trying to incorporate an empty active profile " + profile);

        final ActivityProfile bandPassFiltered = profile.bandPassFilter();
        activeRegions.addAll(bandPassFiltered.createActiveRegions( getActiveRegionExtension(), getMaxRegionSize() ));
        return new ActivityProfile( engine.getGenomeLocParser(), profile.hasPresetRegions() );
    }

    protected final ActivityProfileResult walkerActiveProb(final ActiveRegionWalker<M, T> walker,
                                                           final RefMetaDataTracker tracker, final ReferenceContext refContext,
                                                           final AlignmentContext locus, final GenomeLoc location) {
        if ( walker.hasPresetActiveRegions() ) {
            return new ActivityProfileResult(location, walker.presetActiveRegions.overlaps(location) ? 1.0 : 0.0);
        } else {
            return walker.isActive( tracker, refContext, locus );
        }
    }

    protected ReferenceOrderedView getReferenceOrderedView(final ActiveRegionWalker<M, T> walker,
                                                           final LocusShardDataProvider dataProvider,
                                                           final LocusView locusView) {
        if ( WalkerManager.getWalkerDataSource(walker) != DataSource.REFERENCE_ORDERED_DATA )
            return new ManagingReferenceOrderedView( dataProvider );
        else
            return (RodLocusView)locusView;
    }

    /**
     * Write out each active region to the walker activeRegionOutStream
     *
     * @param walker
     */
    protected void writeActiveRegionsToStream( final ActiveRegionWalker<M, T> walker ) {
        // Just want to output the active regions to a file, not actually process them
        for( final ActiveRegion activeRegion : workQueue ) {
            if( activeRegion.isActive ) {
                walker.activeRegionOutStream.println( activeRegion.getLocation() );
            }
        }
    }
}
