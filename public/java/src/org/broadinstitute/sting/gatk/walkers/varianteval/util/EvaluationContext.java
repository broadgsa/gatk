package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

public final class EvaluationContext {
    // NOTE: must be hashset to avoid O(log n) cost of iteration in the very frequently called apply function
    private final HashSet<VariantEvaluator> evaluationInstances;

    public EvaluationContext(final VariantEvalWalker walker, final Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        evaluationInstances = new HashSet<VariantEvaluator>(evaluationClasses.size());

        for ( final Class<? extends VariantEvaluator> c : evaluationClasses ) {
            try {
                final VariantEvaluator eval = c.newInstance();
                eval.initialize(walker);
                evaluationInstances.add(eval);
            } catch (InstantiationException e) {
                throw new ReviewedStingException("Unable to instantiate eval module '" + c.getSimpleName() + "'", e);
            } catch (IllegalAccessException e) {
                throw new ReviewedStingException("Illegal access error when trying to instantiate eval module '" + c.getSimpleName() + "'", e);
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
                    throw new ReviewedStingException("BUG: Unexpected evaluation order " + evaluation);
            }
        }
    }
}
