package org.broadinstitute.sting.gatk.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 19, 2009
 */
public class NoOriginalQualityScoresFilter implements SamRecordFilter {
    public boolean filterOut(SAMRecord rec) {
        return (rec.getAttribute("OQ") == null);
    }
}
