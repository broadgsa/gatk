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

import edu.mit.broad.sam.SAMRecord;

import java.util.List;

/**
 * Aggregates multiple filters and provides a method for applying them all to a given record with
 * one method call.
 */
public class AggregateFilter implements SamRecordFilter {

    private final List<SamRecordFilter> filters;

    /**
     * Constructor
     * @param filters   the list of filters that this Aggregator applies
     */
    public AggregateFilter(List<SamRecordFilter> filters) {
        this.filters = filters;
    }

    /**
     * Determines whether a SAMRecord matches this filter
     *
     * @param record    the SAMRecord to evaluate
     * @return  true if the SAMRecord matches at least one filter, otherwise false
     */
    public boolean filterOut(SAMRecord record) {
        for (SamRecordFilter filter : filters) {
            if (filter.filterOut(record)) {
                return true;
            }
        }
        return false;
    }
}
