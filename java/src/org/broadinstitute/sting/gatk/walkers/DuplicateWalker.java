package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.util.Arrays;

import net.sf.samtools.SAMRecord;
import net.sf.picard.filter.SamRecordFilter;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires({DataSource.READS,DataSource.REFERENCE})
public abstract class DuplicateWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(GenomeLoc loc, AlignmentContext context, Set<List<SAMRecord>> readSets ) {
        return true;    // We are keeping all the reads
    }

    public abstract MapType map(GenomeLoc loc, AlignmentContext context, Set<List<SAMRecord>> readSets );

    // Given result of map function
    public abstract ReduceType reduceInit();
    public abstract ReduceType reduce(MapType value, ReduceType sum);

    // --------------------------------------------------------------------------------------------------------------
    //
    // read filters
    //
    // --------------------------------------------------------------------------------------------------------------

    public List<SamRecordFilter> getMandatoryReadFilters() {
        SamRecordFilter filter1 = new UnmappedReadFilter();
        SamRecordFilter filter2 = new NotPrimaryAlignmentReadFilter();
        List<SamRecordFilter> x = super.getMandatoryReadFilters();

        x.addAll(Arrays.asList(filter2, filter1));
        return x;

    }
}