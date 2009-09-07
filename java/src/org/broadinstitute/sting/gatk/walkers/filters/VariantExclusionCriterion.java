package org.broadinstitute.sting.gatk.walkers.filters;

import java.util.HashMap;

public interface VariantExclusionCriterion {
    public void initialize(HashMap<String,String> arguments);
    public void compute(VariantContextWindow contextWindow);

    //public boolean isExcludable();
    public double inclusionProbability();

    public String getStudyHeader();

    public String getStudyInfo();
    public boolean useZeroQualityReads();
}
