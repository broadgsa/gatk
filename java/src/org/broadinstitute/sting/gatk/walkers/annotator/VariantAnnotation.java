package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Variation;


public interface VariantAnnotation {

    // return the annotation for the given locus data (return null for no annotation)
    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation);

    // return true if you want to use a context with mapping quality zero reads, false otherwise
    public boolean useZeroQualityReads();

    // return the INFO key
    public String getKeyName();

    // return the description used for the VCF INFO meta field
    public String getDescription();

}
