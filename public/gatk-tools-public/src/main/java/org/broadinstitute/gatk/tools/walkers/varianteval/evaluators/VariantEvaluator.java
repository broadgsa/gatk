/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import htsjdk.variant.variantcontext.VariantContext;

public abstract class VariantEvaluator implements Comparable<VariantEvaluator> {
    private VariantEval walker;
    private final String simpleName;

    protected VariantEvaluator() {
        this.simpleName = getClass().getSimpleName();
    }

    public void initialize(VariantEval walker) {
        this.walker = walker;
    }

    public VariantEval getWalker() {
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
        return eval.getAttributeAsBoolean(VariantEval.IS_SINGLETON_KEY, false);
    }

    public final String getSimpleName() {
        return simpleName;
    }

    @Override
    public int compareTo(final VariantEvaluator variantEvaluator) {
        return getSimpleName().compareTo(variantEvaluator.getSimpleName());
    }

    /**
     * Evaluation modules that override this function to indicate that they support
     * combining the results of two independent collections of eval data into
     * a single meaningful result.  The purpose of this interface is to
     * allow us to cut up the input data into many independent stratifications, and then
     * at the end of the eval run decide which stratifications to combine.  This is
     * important in the case of AC, where you may have thousands of distinct AC
     * values that chop up the number of variants to too small a number of variants,
     * and you'd like to combine the AC values into ranges containing some percent
     * of the data.
     *
     * For example, suppose you have an eval that
     * counts variants in a variable nVariants.  If you want to be able to combine
     * multiple evaluations of this type, overload the combine function
     * with a function that sets this.nVariants += other.nVariants.
     *
     * Add in the appropriate fields of the VariantEvaluator T
     * (of the same type as this object) to the values of this object.
     *
     * The values in this and other are implicitly independent, so that
     * the values can be added together.
     *
     * @param other a VariantEvaluator of the same type of this object
     */
    public void combine(final VariantEvaluator other) {
        throw new ReviewedGATKException(getSimpleName() + " doesn't support combining results, sorry");
    }

    /**
     * Must be overloaded to return true for evaluation modules that support the combine operation
     *
     * @return
     */
    public boolean supportsCombine() {
        return false;
    }
}
