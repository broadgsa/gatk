package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEval;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
