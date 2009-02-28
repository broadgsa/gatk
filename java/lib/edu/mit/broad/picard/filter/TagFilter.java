/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.filter;

import edu.mit.broad.sam.SAMRecord;

import java.util.List;
import java.util.Arrays;

/**
 * Filter class for matching tag attributes in SAMRecords
 */
public class TagFilter implements SamRecordFilter {

    private final String tag;           // The key of the tag to match
    private final List<Object> values;  // The list of matching values

    /**
     * Constructor for a single value
     *
     * @param tag       the key of the tag to match
     * @param value     the value to match
     */
    public TagFilter(String tag, Object value) {
        this.tag = tag;
        this.values = Arrays.asList(value);
    }

    /**
     * Constructor for multiple values
     *
     * @param tag       the key of the tag to match
     * @param values    the matching values
     */
    public TagFilter(String tag, List<Object> values) {
        this.tag = tag;
        this.values = values;
    }

    /**
     * Determines whether a SAMRecord matches this filter
     *
     * @param record    the SAMRecord to evaluate
     * @return  true if the SAMRecord matches the filter, otherwise false
     */
    public boolean filterOut(SAMRecord record) {
        return values.contains(record.getAttribute(tag));
    }
 }
