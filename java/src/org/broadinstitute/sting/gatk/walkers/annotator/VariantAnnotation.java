package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;


public interface VariantAnnotation {

    // return the annotation for the given variation and context split by sample (return null for no annotation)
    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation);

    // return the INFO key
    public String getKeyName();

    // return the description used for the VCF INFO meta field
    public VCFInfoHeaderLine getDescription();

}
