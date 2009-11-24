package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.List;


public class DepthOfCoverage extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        int depth = pileup.getReads().size();
        return new Pair<String, String>("DoC", String.format("%d", depth));
    }

    public boolean useZeroQualityReads() { return false; }
}
