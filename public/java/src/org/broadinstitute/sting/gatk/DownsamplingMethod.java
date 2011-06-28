package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.utils.exceptions.UserException;

/**
 * Describes the method for downsampling reads at a given locus.
 *
 * @author hanna
 * @version 0.1
 */

public class DownsamplingMethod {
    /**
     * Type of downsampling to perform.
     */
    public final DownsampleType type;

    /**
     * Actual downsampling target is specified as an integer number of reads.
     */
    public final Integer toCoverage;

    /**
     * Actual downsampling target is specified as a fraction of total available reads.
     */
    public final Double toFraction;

    /**
     * Expresses no downsampling applied at all.
     */
    public static final DownsamplingMethod NONE = new DownsamplingMethod(DownsampleType.NONE,null,null);

    public DownsamplingMethod(DownsampleType type, Integer toCoverage, Double toFraction) {
        // Do some basic sanity checks on the downsampling parameters passed in.

        // Can't leave toFraction and toCoverage null unless type is experimental naive duplicate eliminator.
        if(type != DownsampleType.NONE && toFraction == null && toCoverage == null)
            throw new UserException.CommandLineException("Must specify either toFraction or toCoverage when downsampling.");

        // Fraction and coverage cannot both be specified.
        if(toFraction != null && toCoverage != null)
            throw new UserException.CommandLineException("Downsampling coverage and fraction are both specified.  Please choose only one.");

        // Experimental by sample downsampling does not work with a fraction of reads.
        if(type == DownsampleType.BY_SAMPLE && toFraction != null)
            throw new UserException.CommandLineException("Cannot downsample to fraction with new EXPERIMENTAL_BY_SAMPLE method");

        this.type = type;
        this.toCoverage = toCoverage;
        this.toFraction = toFraction;
    }
}
