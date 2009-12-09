package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;


public class DepthOfCoverage extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation) {
        int depth = pileup.getReads().size();
        return String.format("%d", depth);
    }

    public String getKeyName() { return VCFRecord.DEPTH_KEY; }

    public String getDescription() { return getKeyName() + ",1,Integer,\"Total Depth (including MQ0 reads)\""; }

    public boolean useZeroQualityReads() { return true; }
}
