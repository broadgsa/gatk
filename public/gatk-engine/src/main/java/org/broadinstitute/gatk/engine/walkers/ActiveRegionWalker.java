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

import com.google.java.contract.Ensures;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.engine.filters.*;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegion;
import org.broadinstitute.gatk.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.gatk.utils.activeregion.ActivityProfileState;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.interval.IntervalSetRule;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Base class for all the Active Region Walkers.
 * User: rpoplin
 * Date: 12/7/11
 */

@By(DataSource.READS)
@Requires({DataSource.READS, DataSource.REFERENCE})
@PartitionBy(PartitionType.READ)
@ActiveRegionTraversalParameters(extension=50,maxRegion=1500)
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, MappingQualityUnavailableFilter.class})
@Downsample(by = DownsampleType.BY_SAMPLE, toCoverage = 1000)
@RemoveProgramRecords
public abstract class ActiveRegionWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    /**
     * If provided, this walker will write out its activity profile (per bp probabilities of being active)
     * to this file in the IGV formatted TAB deliminated output:
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * Intended to make debugging the activity profile calculations easier
     */
    @Output(fullName="activityProfileOut", shortName="APO", doc="Output the raw activity profile results in IGV format", required = false, defaultToStdout = false)
    public PrintStream activityProfileOutStream = null;

    /**
     * If provided, this walker will write out its active and inactive regions
     * to this file in the IGV formatted TAB deliminated output:
     *
     * http://www.broadinstitute.org/software/igv/IGV
     *
     * Intended to make debugging the active region calculations easier
     */
    @Output(fullName="activeRegionOut", shortName="ARO", doc="Output the active region to this IGV formatted file", required = false, defaultToStdout = false)
    public PrintStream activeRegionOutStream = null;

    @Advanced
    @Input(fullName="activeRegionIn", shortName="AR", doc="Use this interval list file as the active regions to process", required = false)
    protected List<IntervalBinding<Feature>> activeRegionBindings = null;

    @Advanced
    @Argument(fullName="activeRegionExtension", shortName="activeRegionExtension", doc="The active region extension; if not provided defaults to Walker annotated default", required = false)
    public Integer activeRegionExtension = null;

    /**
     * For the active region walker to treat all bases as active.  Useful for debugging when you want to force something like
     * the HaplotypeCaller to process a specific interval you provide the GATK
     */
    @Advanced
    @Argument(fullName="forceActive", shortName="forceActive", doc="If provided, all bases will be tagged as active", required = false)
    public boolean forceActive = false;

    @Advanced
    @Argument(fullName="activeRegionMaxSize", shortName="activeRegionMaxSize", doc="The active region maximum size; if not provided defaults to Walker annotated default", required = false)
    public Integer activeRegionMaxSize = null;

    @Advanced
    @Argument(fullName="bandPassSigma", shortName="bandPassSigma", doc="The sigma of the band pass filter Gaussian kernel; if not provided defaults to Walker annotated default", required = false)
    public Double bandPassSigma = null;

    /*
     * For active region limits in ActivityProfile
*   */
    @Hidden
    @Argument(fullName = "maxProbPropagationDistance", shortName = "maxProbPropDist", minValue = 0, doc="Region probability propagation distance beyond it's maximum size.", required = false)
    public Integer maxProbPropagationDistance = 50;

    @Advanced
    @Argument(fullName = "activeProbabilityThreshold", shortName = "ActProbThresh", minValue = 0.0, maxValue = 1.0, doc="Threshold for the probability of a profile state being active.", required = false)
    public Double activeProbThreshold = 0.002;

    private GenomeLocSortedSet presetActiveRegions = null;

    @Override
    public void initialize() {
        if( activeRegionBindings == null ) { return; }
        List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>(0);
        for ( IntervalBinding intervalBinding : activeRegionBindings ) {
            List<GenomeLoc> intervals = intervalBinding.getIntervals(this.getToolkit().getGenomeLocParser());

            if ( intervals.isEmpty() ) {
                logger.warn("The interval file " + intervalBinding.getSource() + " contains no intervals that could be parsed.");
            }

            allIntervals = IntervalUtils.mergeListsBySetOperator(intervals, allIntervals, IntervalSetRule.UNION);
        }

        presetActiveRegions = IntervalUtils.sortAndMergeIntervals(this.getToolkit().getGenomeLocParser(), allIntervals, IntervalMergingRule.ALL);
    }

    /**
     * Does this walker want us to use a set of preset action regions instead of dynamically using the result of isActive?
     * @return true if yes, false if no
     */
    public boolean hasPresetActiveRegions() {
        return presetActiveRegions != null;
    }

    /**
     * Get the set of preset active regions, or null if none were provided
     * @return a set of genome locs specifying fixed active regions requested by the walker, or null if none exist
     */
    public GenomeLocSortedSet getPresetActiveRegions() {
        return presetActiveRegions;
    }

    // Do we actually want to operate on the context?
    public boolean filter(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    public EnumSet<ActiveRegionReadState> desiredReadStates() {
        return EnumSet.of(ActiveRegionReadState.PRIMARY);
    }

    public final boolean wantsNonPrimaryReads() {
        return desiredReadStates().contains(ActiveRegionReadState.NONPRIMARY);
    }

    public boolean wantsExtendedReads() {
        return desiredReadStates().contains(ActiveRegionReadState.EXTENDED);
    }

    public boolean wantsUnmappedReads() {
        return desiredReadStates().contains(ActiveRegionReadState.UNMAPPED);
    }

    // Determine probability of active status over the AlignmentContext
    @Ensures({"result.isActiveProb >= 0.0", "result.isActiveProb <= 1.0"})
    public abstract ActivityProfileState isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context);

    // Map over the ActiveRegion
    public abstract MapType map(final ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker);

    public final GenomeLocSortedSet extendIntervals( final GenomeLocSortedSet intervals, final GenomeLocParser genomeLocParser, IndexedFastaSequenceFile reference ) {
        final int activeRegionExtension = this.getClass().getAnnotation(ActiveRegionTraversalParameters.class).extension();
        final List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>();
        for( final GenomeLoc interval : intervals.toList() ) {
            final int start = Math.max( 1, interval.getStart() - activeRegionExtension );
            final int stop = Math.min( reference.getSequenceDictionary().getSequence(interval.getContig()).getSequenceLength(), interval.getStop() + activeRegionExtension );
            allIntervals.add( genomeLocParser.createGenomeLoc(interval.getContig(), start, stop) );
        }
        return IntervalUtils.sortAndMergeIntervals(genomeLocParser, allIntervals, IntervalMergingRule.ALL);
    }


}
