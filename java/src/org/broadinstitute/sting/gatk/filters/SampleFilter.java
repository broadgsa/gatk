package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.utils.cmdLine.Argument;

public class SampleFilter implements SamRecordFilter {
    @Argument(fullName = "sample_to_keep", shortName = "goodSM", doc="The name of the sample to keep, filtering out all others", required=true)
    private String SAMPLE_TO_KEEP = null;

    public boolean filterOut( final SAMRecord read ) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        return !( readGroup != null && readGroup.getSample().equals( SAMPLE_TO_KEEP ) );
    }
}
