package org.broadinstitute.sting.oneoffprojects.filters;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMRecord;

import java.util.Arrays;
import java.util.HashSet;

/**
 * IF THERE IS NO JAVADOC RIGHT HERE, YELL AT chartl
 *
 * @Author chartl
 * @Date Feb 8, 2010
 */
public class ContaminatedSampleFilter implements SamRecordFilter {

    private final String[] filteredNames = {"NA19562","NA19006","NA19554","NA18985","NA18988","NA19062","NA19559"};

    private HashSet<String> contaminatedSamples = new HashSet<String>( Arrays.asList(filteredNames) );

    public boolean filterOut(SAMRecord read) {
        return contaminatedSamples.contains(read.getReadGroup().getSample());
    }
}
