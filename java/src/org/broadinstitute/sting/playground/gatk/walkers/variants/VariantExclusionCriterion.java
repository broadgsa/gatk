package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;

public interface VariantExclusionCriterion {
    public void initialize(String arguments);

    public void compute(char ref, AlignmentContext context, rodVariants variant);

    //public boolean isExcludable();
    public double inclusionProbability();

    public String getStudyHeader();

    public String getStudyInfo();
    public boolean useZeroQualityReads();
}
