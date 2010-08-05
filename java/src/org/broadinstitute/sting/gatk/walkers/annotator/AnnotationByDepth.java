package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import java.util.Map;



public abstract class AnnotationByDepth implements InfoFieldAnnotation {


    protected int annotationByVariantDepth(final Map<String, Genotype> genotypes, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
        int depth = 0;
        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {

            // we care only about variant calls
            if ( genotype.getValue().isHomRef() )
                continue;

            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context != null )
                depth += context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();
        }

        return depth;
    }


}