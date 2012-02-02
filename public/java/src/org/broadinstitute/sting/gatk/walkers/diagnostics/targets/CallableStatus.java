package org.broadinstitute.sting.gatk.walkers.diagnostics.targets;

/**
 * Short one line description of the walker.
 *
 * @author Mauricio Carneiro
 * @since 2/1/12
 */
public enum CallableStatus {
    /** the reference base was an N, which is not considered callable the GATK */
    REF_N,
    /** the base satisfied the min. depth for calling but had less than maxDepth to avoid having EXCESSIVE_COVERAGE */
    CALLABLE,
    /** absolutely no reads were seen at this locus, regardless of the filtering parameters */
    NO_COVERAGE,
    /** there were less than min. depth bases at the locus, after applying filters */
    LOW_COVERAGE,
    /** more than -maxDepth read at the locus, indicating some sort of mapping problem */
    EXCESSIVE_COVERAGE,
    /** more than --maxFractionOfReadsWithLowMAPQ at the locus, indicating a poor mapping quality of the reads */
    POOR_QUALITY
}
