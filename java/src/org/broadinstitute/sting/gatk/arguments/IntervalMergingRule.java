package org.broadinstitute.sting.gatk.arguments;


/**
 * a class we use to determine the merging rules for intervals passed to the GATK
 */
public enum IntervalMergingRule {
    ALL, // we merge both overlapping intervals and abutting intervals
    OVERLAPPING_ONLY, // We merge intervals that are overlapping, but NOT ones that only abut each other
    NONE; // we merge neither overlapping or abutting intervals, the list of intervals is sorted, but not merged

    public boolean check() {
        if (this.compareTo(NONE) == 0)
            throw new UnsupportedOperationException("We Currently do not support IntervalMergingRule.NONE");
        return true;
    }
}
