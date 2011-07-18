package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

public interface GenotypeAnnotation {

    // return annotations for the given contexts/genotype split by sample
    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext stratifiedContext, VariantContext vc, Genotype g);

    // return the FORMAT keys
    public List<String> getKeyNames();

    // return the descriptions used for the VCF FORMAT meta field
    public List<VCFFormatHeaderLine> getDescriptions();
   
}