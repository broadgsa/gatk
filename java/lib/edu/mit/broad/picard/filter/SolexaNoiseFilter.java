/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.filter;

import edu.mit.broad.picard.util.SequenceUtil;
import edu.mit.broad.sam.SAMRecord;

/**
 * Filter to determine whether a read is "noisy" due to a poly-A run that is a sequencing artifact.
 * Currently we filter out only reads that are composed entirely of As.
 */
public class SolexaNoiseFilter implements SamRecordFilter {

    /**
     * Determines whether a SAMRecord matches this filter
     *
     * @param record    the SAMRecord to evaluate
     * @return  true if the SAMRecord matches the filter, otherwise false
     */
    public boolean filterOut(SAMRecord record) {
        byte sequence[] = record.getReadBases();
        for (byte base : sequence) {
            if (base != 'A' && base != 'a' &&
                !SequenceUtil.isNoCall(base)) {
                return false;
            }
        }
        return true;
    }
}
