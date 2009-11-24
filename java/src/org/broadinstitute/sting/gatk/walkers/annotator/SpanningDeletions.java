package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.List;


public class SpanningDeletions extends StandardVariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {
        int deletions = 0;
        for (Integer offset : pileup.getOffsets() ) {
            if ( offset == -1 )
                deletions++;
        }
        return new Pair<String, String>("SpanningDeletionFraction", String.format("%.2f", (double)deletions/(double)pileup.getReads().size()));
    }

    public boolean useZeroQualityReads() { return false; }
}