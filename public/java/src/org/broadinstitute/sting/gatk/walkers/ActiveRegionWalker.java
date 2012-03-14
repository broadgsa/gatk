package org.broadinstitute.sting.gatk.walkers;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.gatk.filters.FailsVendorQualityCheckFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.activeregion.ActiveRegion;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalSetRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Base class for all the Active Region Walkers.
 * User: rpoplin
 * Date: 12/7/11
 */

@By(DataSource.READS)
@Requires({DataSource.READS, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.READ)
@ActiveRegionExtension(extension=50)
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class})
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

    public boolean wantsNonPrimaryReads() {
        return false;
    }

    // Determine probability of active status over the AlignmentContext
    public abstract double isActive(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context);

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
