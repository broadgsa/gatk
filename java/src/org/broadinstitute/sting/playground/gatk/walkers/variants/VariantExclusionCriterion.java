package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;

public interface VariantExclusionCriterion {
    public void initialize(String arguments);

    public boolean exclude(char ref, LocusContext context, rodVariants variant);
    public boolean useZeroQualityReads();
}
