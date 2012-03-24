package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public abstract class VariantEvaluator {
    private VariantEvalWalker walker;

    public void initialize(VariantEvalWalker walker) {
        this.walker = walker;
    }

    public VariantEvalWalker getWalker() {
        return walker;
    }

    public abstract boolean enabled();

    // Should return the number of VariantContexts expected as inputs to update.  Can be 1 or 2
    public abstract int getComparisonOrder();

    // called at all sites, regardless of eval context itself; useful for counting processed bases
    // No longer available.  The processed bp is kept in VEW itself for performance reasons
    //    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

    public String update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    public void finalizeEvaluation() {}

    protected double rate(long n, long d) {
        return n / (1.0 * Math.max(d, 1));
    }

    protected long inverseRate(long n, long d) {
        return n == 0 ? 0 : d / Math.max(n, 1);
    }

    protected double ratio(long num, long denom) {
        return ((double)num) / (Math.max(denom, 1));
    }

    /**
     * Returns true if the variant in vc was a singleton in the original input evaluation
     * set, regardless of variant context subsetting that has occurred.
     * @param eval the VariantContext being assessed for this previous status as a singleton
     * @return true if eval was originally a singleton site
     */
    protected static boolean variantWasSingleton(final VariantContext eval) {
        return eval.getAttributeAsBoolean(VariantEvalWalker.IS_SINGLETON_KEY, false);
    }

    /**
     * Convenience function that formats the novelty rate as a %.2f string
     *
     * @param known number of variants from all that are known
     * @param all number of all variants
     * @return a String novelty rate, or NA if all == 0
     */
    protected static String formattedNoveltyRate(final int known, final int all) {
        return formattedPercent(all - known, all);
    }

    /**
     * Convenience function that formats the novelty rate as a %.2f string
     *
     * @param x number of objects part of total that meet some criteria
     * @param total count of all objects, including x
     * @return a String percent rate, or NA if total == 0
     */
    protected static String formattedPercent(final int x, final int total) {
        return total == 0 ? "NA" : String.format("%.2f", x / (1.0*total));
    }

    /**
     * Convenience function that formats a ratio as a %.2f string
     *
     * @param num  number of observations in the numerator
     * @param denom number of observations in the denumerator
     * @return a String formatted ratio, or NA if all == 0
     */
    protected static String formattedRatio(final int num, final int denom) {
        return denom == 0 ? "NA" : String.format("%.2f", num / (1.0 * denom));
    }
}
