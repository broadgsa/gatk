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

import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;

public class NewEvaluationContext extends HashMap<VariantStratifier, String> {
    public TreeMap<String, VariantEvaluator> evaluationInstances;

    public String toString() {
        String value = "";

        for ( VariantStratifier key : this.keySet() ) {
            value += "\t" + key.getClass().getSimpleName() + ":" + this.get(key) + "\n";
        }

        return value;
    }

    public void addEvaluationClassList(VariantEvalWalker walker, StateKey stateKey, Set<Class<? extends VariantEvaluator>> evaluationClasses) {
        evaluationInstances = new TreeMap<String, VariantEvaluator>();

        for ( Class<? extends VariantEvaluator> c : evaluationClasses ) {
            try {
                VariantEvaluator eval = c.newInstance();
                eval.initialize(walker);

                if (eval.stateIsApplicable(stateKey)) {
                    evaluationInstances.put(c.getSimpleName(), eval);
                }
            } catch (InstantiationException e) {
                throw new StingException("Unable to instantiate eval module '" + c.getSimpleName() + "'");
            } catch (IllegalAccessException e) {
                throw new StingException("Illegal access error when trying to instantiate eval module '" + c.getSimpleName() + "'");
            }
        }
    }

    public TreeMap<String, VariantEvaluator> getEvaluationClassList() {
        return evaluationInstances;
    }

    public void apply(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantContext comp, VariantContext eval) {
        for ( VariantEvaluator evaluation : evaluationInstances.values() ) {
            // we always call update0 in case the evaluation tracks things like number of bases covered

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
                    //if (eval != null) {
                        evaluation.update2(eval, comp, tracker, ref, context);
                    //}

                    break;
                default:
                    throw new ReviewedStingException("BUG: Unexpected evaluation order " + evaluation);
            }
        }
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( VariantEvaluator evaluation : evaluationInstances.values() ) {
            evaluation.update0(tracker, ref, context);
        }
    }
}
