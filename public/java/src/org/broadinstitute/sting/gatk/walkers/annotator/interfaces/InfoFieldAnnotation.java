package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

public abstract class InfoFieldAnnotation extends VariantAnnotatorAnnotation {
    // return annotations for the given contexts split by sample
    public abstract Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker,
                                                 ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc);

    // return the descriptions used for the VCF INFO meta field
    public abstract List<VCFInfoHeaderLine> getDescriptions();
}