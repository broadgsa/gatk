package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.genotype.vcf.VCFFormatHeaderLine;

import java.util.Map;

public interface GenotypeAnnotation {

    // annotate the given record for the given variation and context split by sample
    public void annotateContext(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc);

    // return the FORMAT key
    public String getKeyName();

    // return the description used for the VCF FORMAT meta field
    public VCFFormatHeaderLine getDescription();
   
}