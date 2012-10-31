package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

// TODO -- make this an abstract class when we move away from InfoFieldAnnotation
public interface ActiveRegionBasedAnnotation extends AnnotationType {
    // return annotations for the given contexts split by sample and then read likelihood
    public abstract Map<String, Object> annotate(final Map<String,PerReadAlleleLikelihoodMap> stratifiedContexts, final VariantContext vc);

    // return the descriptions used for the VCF INFO meta field
    public abstract List<VCFInfoHeaderLine> getDescriptions();
}