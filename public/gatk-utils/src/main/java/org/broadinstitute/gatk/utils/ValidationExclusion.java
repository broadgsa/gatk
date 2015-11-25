/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils;

import org.broadinstitute.gatk.utils.commandline.EnumerationArgumentDefault;

import java.util.ArrayList;
import java.util.List;


public class ValidationExclusion {
    // our validation options

    public enum TYPE {
        ALLOW_N_CIGAR_READS,    // ignore the presence of N operators in CIGARs: do not blow up and process reads that contain one or more N operators.
                                // This exclusion does not have effect on reads that get filtered {@see MalformedReadFilter}.
        ALLOW_UNINDEXED_BAM,        // allow bam files that do not have an index; we'll traverse them using monolithic shard
        ALLOW_UNSET_BAM_SORT_ORDER, // assume that the bam is sorted, even if the SO (sort-order) flag is not set
        NO_READ_ORDER_VERIFICATION, // do not validate that the reads are in order as we take them from the bam file
        ALLOW_SEQ_DICT_INCOMPATIBILITY, // allow dangerous, but not fatal, sequence dictionary incompatibilities
        LENIENT_VCF_PROCESSING,         // allow non-standard values for standard VCF header lines.  Don't worry about size differences between header and values, etc.
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

    public static boolean lenientVCFProcessing(final TYPE val) {
        return val == TYPE.ALL
                || val == TYPE.LENIENT_VCF_PROCESSING;
    }
}
