package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.traversals.TraversalStatistics;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentReadFilter;
import org.broadinstitute.sting.gatk.filters.DuplicateReadFilter;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;

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
public abstract class LocusWalker<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);

    // --------------------------------------------------------------------------------------------------------------
    //
    // mandatory read filters
    //
    // --------------------------------------------------------------------------------------------------------------

    public List<SamRecordFilter> getMandatoryReadFilters() {
//        if ( false ) {
//            SamRecordFilter filter = new LocusStreamFilterFunc();
//            return Arrays.asList(filter);
//        } else {
        SamRecordFilter filter1 = new UnmappedReadFilter();
        SamRecordFilter filter2 = new NotPrimaryAlignmentReadFilter();
        SamRecordFilter filter3 = new DuplicateReadFilter();
        List<SamRecordFilter> x = super.getMandatoryReadFilters();
        x.addAll(Arrays.asList(filter3, filter2, filter1));
//        }
        return x;
    }

    /**
     * Class to filter out un-handle-able reads from the stream.  We currently are skipping
     * unmapped reads, non-primary reads, unaligned reads, and duplicate reads.
     */
    private static class LocusStreamFilterFunc implements SamRecordFilter {
        SAMRecord lastRead = null;
        public boolean filterOut(SAMRecord rec) {
            boolean result = false;
            if (rec.getReadUnmappedFlag()) {
                TraversalStatistics.nUnmappedReads++;
                result = true;
                //why = "Unmapped";
            } else if (rec.getNotPrimaryAlignmentFlag()) {
                TraversalStatistics.nNotPrimary++;
                result = true;
//                why = "Not Primary";
            } else if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                TraversalStatistics.nBadAlignments++;
                result = true;
//                why = "No alignment start";
            } else if (rec.getDuplicateReadFlag()) {
                TraversalStatistics.nDuplicates++;
                result = true;
//                why = "Duplicate reads";
            }
            else {
                result = false;
            }

            if (result) {
                TraversalStatistics.nSkippedReads++;
                //System.out.printf("  [filter] %s => %b %s", rec.getReadName(), result, why);
            } else {
                TraversalStatistics.nReads++;
            }

            return result;
        }
    }
}
