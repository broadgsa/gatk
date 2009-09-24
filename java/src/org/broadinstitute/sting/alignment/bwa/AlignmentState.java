package org.broadinstitute.sting.alignment.bwa;

/**
 * The current state of an alignment.
 *
 * @author mhanna
 * @version 0.1
 */
public enum AlignmentState {
    MATCH_MISMATCH,
    INSERTION,
    DELETION
}
