package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

public abstract class VariantEvaluator implements Comparable<VariantEvaluator> {
    private VariantEvalWalker walker;
    private final String simpleName;

    protected VariantEvaluator() {
        this.simpleName = getClass().getSimpleName();
    }

    public void initialize(VariantEvalWalker walker) {
        this.walker = walker;
    }

    public VariantEvalWalker getWalker() {
        return walker;
    }

    // Should return the number of VariantContexts expected as inputs to update.  Can be 1 or 2
    public abstract int getComparisonOrder();

    // called at all sites, regardless of eval context itself; useful for counting processed bases
    // No longer available.  The processed bp is kept in VEW itself for performance reasons
    //    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

    public void update1(VariantContext eval, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
    }

    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
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

    public String getSimpleName() {
        return simpleName;
    }

    @Override
    public int compareTo(final VariantEvaluator variantEvaluator) {
        return getSimpleName().compareTo(variantEvaluator.getSimpleName());
    }
}
