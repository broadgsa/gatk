package org.broadinstitute.sting.gatk.datasources.simpleDataSources;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.List;

/**
 * Filters reads out of a data stream that don't overlap with the given list of locations.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalOverlappingFilter implements SamRecordFilter {
    /**
     * The list of locations containing reads to keep.
     */
    private final List<GenomeLoc> intervals;

    public IntervalOverlappingFilter(List<GenomeLoc> intervals) {
        this.intervals = intervals;
    }

    /**
     * Filter out this record if it doesn't appear in the interval list.
     * @param read The read to examine.
     * @return True to filter the read out.  False otherwise.
     */
    public boolean filterOut(SAMRecord read) {
        GenomeLoc readLocation = GenomeLocParser.createGenomeLoc(read);
        for(GenomeLoc interval: intervals) {
            if(interval.overlapsP(readLocation))
                return false;
        }
        return true;
    }

}
