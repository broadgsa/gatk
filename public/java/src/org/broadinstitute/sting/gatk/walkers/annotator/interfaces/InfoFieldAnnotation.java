package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

public abstract class InfoFieldAnnotation extends VariantAnnotatorAnnotation {
    // return annotations for the given contexts split by sample
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc) {
        return annotate(tracker, walker, ref, stratifiedContexts, vc, null);
    }

    public Map<String, Object> annotate(Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap, VariantContext vc) {
        return annotate(null, null, null, null, vc, perReadAlleleLikelihoodMap);
    }


    public abstract Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                                 final AnnotatorCompatible walker,
                                                 final ReferenceContext ref,
                                                 final Map<String, AlignmentContext> stratifiedContexts,
                                                 final VariantContext vc,
                                                 final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap);

    // return the descriptions used for the VCF INFO meta field
    public abstract List<VCFInfoHeaderLine> getDescriptions();
}