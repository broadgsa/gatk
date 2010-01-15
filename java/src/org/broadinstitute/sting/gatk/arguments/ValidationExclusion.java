package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.utils.cmdLine.EnumerationArgumentDefault;

import java.util.ArrayList;
import java.util.List;


/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 * 
 * @author aaron 
 * 
 * Class ValidationExclusion
 *
 * a class for containing the exclusions from validation that the user
 * wants.
 */
public class ValidationExclusion {

    // our validation options

    public enum TYPE {
        ALLOW_UNINDEXED_BAM,        // allow bam files that do not have an index; we'll traverse them using monolithic shard
        ALLOW_EMPTY_INTERVAL_LIST,  // allow the user to pass in an empty interval list
        ALLOW_UNSET_BAM_SORT_ORDER, // assume that the bam is sorted, even if the SO (sort-order) flag is not set
        NO_READ_ORDER_VERIFICATION, // do not validate that the reads are in order as we take them from the bam file
        @EnumerationArgumentDefault // set the ALL value to the default value, so if they specify just -U, we get the ALL
        ALL                         // do not check for all of the above conditions, DEFAULT
    }

    // a storage for the passed in exclusions
    List<TYPE> exclusions = new ArrayList<TYPE>();

    public ValidationExclusion(List<TYPE> exclusionsList) {
        exclusions.addAll(exclusionsList);
    }

    public ValidationExclusion() {}
    
    /**
     * do we contain the exclusion specified, or were we set to ALL
     * @param t the exclusion case to test for
     * @return true if we contain the exclusion or if we're set to ALL, false otherwise
     */
    public boolean contains(TYPE t) {
        return (exclusions.contains(TYPE.ALL) || exclusions.contains(t));
    }
}
