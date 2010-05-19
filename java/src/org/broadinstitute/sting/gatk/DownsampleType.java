package org.broadinstitute.sting.gatk;

/**
 * Type of downsampling method to invoke.
 *
 * @author hanna
 * @version 0.1
 */

public enum DownsampleType {
    NONE,
    ALL_READS,
    EXPERIMENTAL_BY_SAMPLE,
    EXPERIMENTAL_NAIVE_DUPLICATE_ELIMINATOR
}
