package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS,DataSource.REFERENCE})
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentFilter.class})
public abstract class DuplicateWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(GenomeLoc loc, AlignmentContext context, Set<List<GATKSAMRecord>> readSets ) {
        return true;    // We are keeping all the reads
    }

    public abstract MapType map(GenomeLoc loc, AlignmentContext context, Set<List<GATKSAMRecord>> readSets );

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);
}