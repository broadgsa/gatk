package org.broadinstitute.sting.utils.interval;

/**
 * set operators for combining lists of intervals
 */
public enum IntervalSetRule {
    /** Take the union of all intervals */
    UNION,
    /** Take the intersection of intervals (the subset that overlaps all intervals specified) */
    INTERSECTION;
}
