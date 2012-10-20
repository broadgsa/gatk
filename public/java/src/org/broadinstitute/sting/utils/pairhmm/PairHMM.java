package org.broadinstitute.sting.utils.pairhmm;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

/**
 * Created with IntelliJ IDEA.
 * User: rpoplin
 * Date: 10/16/12
 */

public abstract class PairHMM {
    protected static final Byte MAX_CACHED_QUAL = Byte.MAX_VALUE;
    protected static final byte DEFAULT_GOP = (byte) 45;
    protected static final byte DEFAULT_GCP = (byte) 10;

    public enum HMM_IMPLEMENTATION {
        /* Very slow implementation which uses very accurate log10 sum functions. Only meant to be used as a reference test implementation */
        EXACT,
        /* PairHMM as implemented for the UnifiedGenotyper. Uses log10 sum functions accurate to only 1E-4 */
        ORIGINAL,
        /* Optimized version of the PairHMM which caches per-read computations */
        CACHING,
        /* Optimized version of the PairHMM which caches per-read computations and operations in real space to avoid costly sums of log10'ed likelihoods */
        LOGLESS_CACHING
    }

    protected double[][] matchMetricArray = null;
    protected double[][] XMetricArray = null;
    protected double[][] YMetricArray = null;

    public abstract void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH );

    @Requires({"readBases.length == readQuals.length", "readBases.length == insertionGOP.length", "readBases.length == deletionGOP.length",
               "readBases.length == overallGCP.length", "matchMetricArray!=null", "XMetricArray!=null", "YMetricArray!=null"})
    @Ensures({"!Double.isInfinite(result)", "!Double.isNaN(result)"}) // Result should be a proper log10 likelihood
    public abstract double computeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                     final byte[] readBases,
                                                                     final byte[] readQuals,
                                                                     final byte[] insertionGOP,
                                                                     final byte[] deletionGOP,
                                                                     final byte[] overallGCP,
                                                                     final int hapStartIndex,
                                                                     final boolean recacheReadValues );
}
