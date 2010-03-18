package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFFormatHeaderLine;

import java.util.Map;

public interface GenotypeAnnotation {

    // return annotations for the given contexts/genotype split by sample
    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, StratifiedAlignmentContext stratifiedContext, VariantContext vc, Genotype g);

    // return the FORMAT key
    public String getKeyName();

    // return the description used for the VCF FORMAT meta field
    public VCFFormatHeaderLine getDescription();
   
}