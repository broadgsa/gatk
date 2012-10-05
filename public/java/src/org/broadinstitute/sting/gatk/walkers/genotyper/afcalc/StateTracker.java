package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

/**
 * Keeps track of the best state seen by the exact model and the max states to visit
 * allowing us to abort the search before we visit the entire matrix of AC x samples
 */
final class StateTracker {
    private final static double MAX_LOG10_ERROR_TO_STOP_EARLY = 6; // we want the calculation to be accurate to 1 / 10^6

    final private int[] maxACsToConsider;

    private ExactACcounts ACsAtMax = null;
    private double maxLog10L = Double.NEGATIVE_INFINITY;

    public StateTracker() {
        this(null);
    }

    public StateTracker(final int[] maxACsToConsider) {
        this.maxACsToConsider = maxACsToConsider;
    }

    /**
     * Update the maximum log10L seen, if log10LofKs is higher, and the corresponding ACs of this state
     *
     * @param log10LofKs the likelihood of our current configuration state
     */
    public void update(final double log10LofKs, final ExactACcounts ACs) {
        if ( log10LofKs > getMaxLog10L()) {
            this.setMaxLog10L(log10LofKs);
            this.ACsAtMax = ACs;
        }
    }

    /**
     * Is the likelihood of configuration K too low to consider, related to the
     * maximum likelihood seen already?
     *
     * @param log10LofK the log10 likelihood of the configuration we're considering analyzing
     * @return true if the configuration cannot meaningfully contribute to our likelihood sum
     */
    public boolean tooLowLikelihood(final double log10LofK) {
        return log10LofK < getMaxLog10L() - MAX_LOG10_ERROR_TO_STOP_EARLY;
    }

    /**
     * Are all ACs in otherACs less than or equal to their corresponding ACs in the maxACsToConsider?
     *
     * @param otherACs the set of otherACs that we want to know if we should consider analyzing
     * @return true if otherACs is a state worth considering, or false otherwise
     */
    public boolean withinMaxACs(final ExactACcounts otherACs) {
        if ( maxACsToConsider == null )
            return true;

        final int[] otherACcounts = otherACs.getCounts();

        for ( int i = 0; i < maxACsToConsider.length; i++ ) {
            // consider one more than the max AC to collect a bit more likelihood mass
            if ( otherACcounts[i] > maxACsToConsider[i] + 1 )
                return false;
        }

        return true;
    }

    /**
     * returns true iff all ACs in this object are less than or equal to their corresponding ACs in the provided set
     */
    public boolean isLowerAC(final ExactACcounts otherACs) {
        if ( ACsAtMax == null )
            return true;

        final int[] myACcounts = this.ACsAtMax.getCounts();
        final int[] otherACcounts = otherACs.getCounts();

        for ( int i = 0; i < myACcounts.length; i++ ) {
            if ( myACcounts[i] > otherACcounts[i] )
                return false;
        }

        return true;
    }

    public boolean abort( final double log10LofK, final ExactACcounts ACs ) {
        return tooLowLikelihood(log10LofK) && isLowerAC(ACs);
    }

    public double getMaxLog10L() {
        return maxLog10L;
    }

    public void setMaxLog10L(double maxLog10L) {
        this.maxLog10L = maxLog10L;
    }
}
