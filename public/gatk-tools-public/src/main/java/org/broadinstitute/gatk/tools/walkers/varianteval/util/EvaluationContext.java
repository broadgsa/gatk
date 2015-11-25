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

package org.broadinstitute.gatk.tools.walkers.varianteval.util;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.gatk.tools.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public final class EvaluationContext {
    // NOTE: must be hashset to avoid O(log n) cost of iteration in the very frequently called apply function
    final VariantEval walker;
    private final ArrayList<VariantEvaluator> evaluationInstances;
    private final Set<Class<? extends VariantEvaluator>> evaluationClasses;

    public EvaluationContext(final VariantEval walker, final Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        this(walker, evaluationClasses, true);
    }

    private EvaluationContext(final VariantEval walker, final Set<Class<? extends VariantEvaluator>> evaluationClasses, final boolean doInitialize) {
        this.walker = walker;
        this.evaluationClasses = evaluationClasses;
        this.evaluationInstances = new ArrayList<VariantEvaluator>(evaluationClasses.size());

        for ( final Class<? extends VariantEvaluator> c : evaluationClasses ) {
            try {
                final VariantEvaluator eval = c.newInstance();
                if ( doInitialize ) eval.initialize(walker);
                evaluationInstances.add(eval);
            } catch (InstantiationException e) {
                throw new ReviewedGATKException("Unable to instantiate eval module '" + c.getSimpleName() + "'", e);
            } catch (IllegalAccessException e) {
                throw new ReviewedGATKException("Illegal access error when trying to instantiate eval module '" + c.getSimpleName() + "'", e);
            }
        }
    }

    /**
     * Returns a sorted set of VariantEvaluators
     *
     * @return
     */
    public final TreeSet<VariantEvaluator> getVariantEvaluators() {
        return new TreeSet<VariantEvaluator>(evaluationInstances);
    }

    public final void apply(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantContext comp, VariantContext eval) {
        for ( final VariantEvaluator evaluation : evaluationInstances ) {
            // the other updateN methods don't see a null context
            if ( tracker == null )
                continue;

            // now call the single or paired update function
            switch ( evaluation.getComparisonOrder() ) {
                case 1:
                    if (eval != null) {
                        evaluation.update1(eval, tracker, ref, context);
                    }
                    break;
                case 2:
                    evaluation.update2(eval, comp, tracker, ref, context);
                    break;
                default:
                    throw new ReviewedGATKException("BUG: Unexpected evaluation order " + evaluation);
            }
        }
    }

    public void combine(final EvaluationContext rhs) {
        for ( int i = 0; i < evaluationInstances.size(); i++ )
            evaluationInstances.get(i).combine(rhs.evaluationInstances.get(i));
    }

    public final static EvaluationContextCombiner COMBINER = new EvaluationContext.EvaluationContextCombiner();
    private static class EvaluationContextCombiner implements StratificationManager.Combiner<EvaluationContext> {
        @Override
        public EvaluationContext combine(EvaluationContext lhs, final EvaluationContext rhs) {
            if ( lhs == null )
                lhs = new EvaluationContext(rhs.walker, rhs.evaluationClasses, false);
            lhs.combine(rhs);
            return lhs;
        }
    }
}
