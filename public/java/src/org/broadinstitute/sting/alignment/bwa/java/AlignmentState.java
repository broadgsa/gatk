package org.broadinstitute.sting.alignment.bwa.java;

/**
 * The state of a given base in the alignment.
 *
 * @author mhanna
 * @version 0.1
 */
public enum AlignmentState {
    MATCH_MISMATCH,
    INSERTION,
    DELETION
}
