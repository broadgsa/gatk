package org.broadinstitute.sting.playground.gatk.walkers.variants;

public interface VariantExclusionCriterion {
    public void initialize(String arguments);

    public void compute(VariantContextWindow contextWindow);

    //public boolean isExcludable();
    public double inclusionProbability();

    public String getStudyHeader();

    public String getStudyInfo();
    public boolean useZeroQualityReads();
}
