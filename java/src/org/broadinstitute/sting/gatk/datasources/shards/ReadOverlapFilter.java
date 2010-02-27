package org.broadinstitute.sting.gatk.datasources.shards;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;

/**
 * Filters reads out of a data stream that don't overlap with the given list of locations.
 *
 * @author mhanna
 * @version 0.1
 */
public class ReadOverlapFilter implements SamRecordFilter {
    /**
     * The list of locations containing reads to keep.
     */
    private final List<GenomeLoc> intervals;

    public ReadOverlapFilter(List<GenomeLoc> intervals) {
        this.intervals = intervals;
    }

    /**
     * Filter out this record if it doesn't appear in the interval list.
     * @param read The read to examine.
     * @return True to filter the read out.  False otherwise.
     */
    public boolean filterOut(SAMRecord read) {
        for(GenomeLoc interval: intervals) {
            if((read.getAlignmentStart() >= interval.getStart() && read.getAlignmentStart() <= interval.getStop()) ||
               (read.getAlignmentEnd() >= interval.getStart() && read.getAlignmentEnd() <= interval.getStop()) ||
               (read.getAlignmentStart() < interval.getStart() && read.getAlignmentEnd() > interval.getStop()))
                return false;
        }
        return true;
    }

}
