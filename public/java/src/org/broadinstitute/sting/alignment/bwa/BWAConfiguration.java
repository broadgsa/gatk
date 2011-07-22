package org.broadinstitute.sting.alignment.bwa;

/**
 * Configuration for the BWA/C aligner.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWAConfiguration {
    /**
     * The maximum edit distance used by BWA.
     */
    public Float maximumEditDistance = null;

    /**
     * How many gap opens are acceptable within this alignment?
     */
    public Integer maximumGapOpens = null;

    /**
     * How many gap extensions are acceptable within this alignment?
     */
    public Integer maximumGapExtensions = null;

    /**
     * Do we disallow indels within a certain range from the start / end?
     */
    public Integer disallowIndelWithinRange = null;

    /**
     * What is the scoring penalty for a mismatch?
     */
    public Integer mismatchPenalty = null;

    /**
     * What is the scoring penalty for a gap open?
     */
    public Integer gapOpenPenalty = null;

    /**
     * What is the scoring penalty for a gap extension?
     */
    public Integer gapExtensionPenalty = null;
}
