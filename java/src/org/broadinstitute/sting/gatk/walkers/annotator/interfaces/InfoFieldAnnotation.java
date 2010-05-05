package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

import java.util.Map;

public interface InfoFieldAnnotation {

    // return annotations for the given contexts split by sample
    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc);

    // return the INFO key
    public String getKeyName();

    // return the description used for the VCF INFO meta field
    public VCFInfoHeaderLine getDescription();

}