package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.List;


public class DepthOfCoverage extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        int depth = pileup.getReads().size();
        return new Pair<String, String>(getKeyName(), String.format("%d", depth));
    }

    public String getKeyName() { return "DP"; }

    public String getDescription() { return "DP,1,Integer,\"Total Depth\""; }

    public boolean useZeroQualityReads() { return false; }
}
