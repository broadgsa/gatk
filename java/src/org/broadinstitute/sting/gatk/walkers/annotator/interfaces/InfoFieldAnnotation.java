package org.broadinstitute.sting.gatk.walkers.annotator.interfaces;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;

import java.util.Map;
import java.util.List;

public interface InfoFieldAnnotation {

    // return annotations for the given contexts split by sample
    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc);

    // return the INFO keys
    public List<String> getKeyNames();

    // return the descriptions used for the VCF INFO meta field
    public List<VCFInfoHeaderLine> getDescriptions();

}