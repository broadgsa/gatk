package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.utils.StingException;

/**
 * Configuration for the BWA/C aligner.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWACConfiguration {
    /**
     * ANN (?) file name.
     */
    public String annFileName = null;

    /**
     * AMB (?) file name.
     */
    public String ambFileName = null;

    /**
     * Packed reference sequence file.
     */
    public String pacFileName = null;

    /**
     * Forward BWT file.
     */
    public String forwardBWTFileName = null;

    /**
     * Forward suffix array file.
     */
    public String forwardSAFileName = null;

    /**
     * Reverse BWT file.
     */
    public String reverseBWTFileName = null;

    /**
     * Reverse suffix array file.
     */
    public String reverseSAFileName = null;
    
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

    /**
     * Create a new BWA configuration file using the given prefix.
     * @param prefix Prefix to use when creating the configuration.  Must not be null.
     */
    public BWACConfiguration(String prefix) {
        if(prefix == null)
            throw new StingException("Prefix must not be null.");
        annFileName = prefix + ".ann";
        ambFileName = prefix + ".amb";
        pacFileName = prefix + ".pac";
        forwardBWTFileName = prefix + ".bwt";
        forwardSAFileName = prefix + ".sa";
        reverseBWTFileName = prefix + ".rbwt";
        reverseSAFileName = prefix + ".rsa";
    }

}
