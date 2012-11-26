package org.broadinstitute.sting.gatk.walkers;

import com.google.java.contract.Ensures;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.activeregion.ActiveRegionReadState;
import org.broadinstitute.sting.utils.activeregion.ActivityProfileResult;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

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
@ActiveRegionExtension(extension=50,maxRegion=1500)
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, MappingQualityUnavailableFilter.class})
@RemoveProgramRecords
public abstract class ActiveRegionWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {

    @Output(fullName="activeRegionOut", shortName="ARO", doc="Output the active region to this interval list file", required = false)
    public PrintStream activeRegionOutStream = null;

    @Input(fullName="activeRegionIn", shortName="AR", doc="Use this interval list file as the active regions to process", required = false)
    protected List<IntervalBinding<Feature>> activeRegionBindings = null;

    public GenomeLocSortedSet presetActiveRegions = null;

    public boolean hasPresetActiveRegions() {
        return presetActiveRegions != null;
    }

    @Override
    public void initialize() {
        if( activeRegionBindings == null ) { return; }
        List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>(0);
        for ( IntervalBinding intervalBinding : activeRegionBindings ) {
            List<GenomeLoc> intervals = intervalBinding.getIntervals(this.getToolkit());

            if ( intervals.isEmpty() ) {
                logger.warn("The interval file " + intervalBinding.getSource() + " contains no intervals that could be parsed.");
            }

            allIntervals = IntervalUtils.mergeListsBySetOperator(intervals, allIntervals, IntervalSetRule.UNION);
        }

        presetActiveRegions = IntervalUtils.sortAndMergeIntervals(this.getToolkit().getGenomeLocParser(), allIntervals, IntervalMergingRule.ALL);
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
    public abstract ActivityProfileResult isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context);

    // Map over the ActiveRegion
    public abstract MapType map(final ActiveRegion activeRegion, final RefMetaDataTracker metaDataTracker);

    public final GenomeLocSortedSet extendIntervals( final GenomeLocSortedSet intervals, final GenomeLocParser genomeLocParser, IndexedFastaSequenceFile reference ) {
        final int activeRegionExtension = this.getClass().getAnnotation(ActiveRegionExtension.class).extension();
        final List<GenomeLoc> allIntervals = new ArrayList<GenomeLoc>();
        for( final GenomeLoc interval : intervals.toList() ) {
            final int start = Math.max( 1, interval.getStart() - activeRegionExtension );
            final int stop = Math.min( reference.getSequenceDictionary().getSequence(interval.getContig()).getSequenceLength(), interval.getStop() + activeRegionExtension );
            allIntervals.add( genomeLocParser.createGenomeLoc(interval.getContig(), start, stop) );
        }
        return IntervalUtils.sortAndMergeIntervals(genomeLocParser, allIntervals, IntervalMergingRule.ALL);
    }
}
