package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.List;


public class SpanningDeletions extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        int deletions = pileup.getNumberOfDeletions();
        return new Pair<String, String>(getKeyName(), String.format("%.2f", (double)deletions/(double)pileup.size()));
    }

    public String getKeyName() { return "Dels"; }

    public String getDescription() { return "Dels,1,Float,\"Fraction of Reads Containing Spanning Deletions\""; }

    public boolean useZeroQualityReads() { return false; }
}