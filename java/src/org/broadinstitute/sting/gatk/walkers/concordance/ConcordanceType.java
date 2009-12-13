package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;
import java.util.Set;

public interface ConcordanceType {

    public void initialize(Map<String,String> args, Set<String> samples);
    public String computeConcordance(Map<String, Genotype> samplesToRecords, ReferenceContext ref);
    public String getInfoName();
    public VCFInfoHeaderLine getInfoDescription();
}