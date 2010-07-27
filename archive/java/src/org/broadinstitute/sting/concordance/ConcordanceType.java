package org.broadinstitute.sting.gatk.walkers.concordance;

import org.broad.tribble.vcf.VCFGenotypeRecord;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.Map;
import java.util.Set;

public interface ConcordanceType {

    public void initialize(Map<String,String> args, Set<String> samples);
    public String computeConcordance(Map<String, VCFGenotypeRecord> samplesToRecords, ReferenceContext ref);
    public String getInfoName();
    public VCFInfoHeaderLine getInfoDescription();
}