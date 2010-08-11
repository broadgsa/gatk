package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorFilter;
import net.sf.picard.filter.SamRecordFilter;

import java.util.List;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.READS)
@Requires({DataSource.READS,DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentReadFilter.class,DuplicateReadFilter.class,FailsVendorQualityCheckReadFilter.class})
public abstract class LocusWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);

    /**
     * Returns the set of locus iterator discards that this walker wants the engine to discard automatically
     *
     * By default, locus walkers ignore adaptor bases but still see the both bases in the overlapping but non-adaptor
     * parts of the reads.
     * @return
     */
    public List<LocusIteratorFilter> getDiscards() {
        LocusIteratorFilter filter = new InAdaptorFilter();
        return Arrays.asList(filter);
    }
}
