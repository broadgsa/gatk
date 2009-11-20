package org.broadinstitute.sting.alignment.bwa;

import org.broadinstitute.sting.utils.StingException;

/**
 * Support files for BWT.
 *
 * @author mhanna
 * @version 0.1
 */
public class BWTFiles {
    /**
     * ANN (?) file name.
     */
    public final String annFileName;

    /**
     * AMB (?) file name.
     */
    public final String ambFileName;

    /**
     * Packed reference sequence file.
     */
    public final String pacFileName;

    /**
     * Forward BWT file.
     */
    public final String forwardBWTFileName;

    /**
     * Forward suffix array file.
     */
    public final String forwardSAFileName;

    /**
     * Reverse BWT file.
     */
    public final String reverseBWTFileName;

    /**
     * Reverse suffix array file.
     */
    public final String reverseSAFileName;

    /**
     * Create a new BWA configuration file using the given prefix.
     * @param prefix Prefix to use when creating the configuration.  Must not be null.
     */
    public BWTFiles(String prefix) {
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
