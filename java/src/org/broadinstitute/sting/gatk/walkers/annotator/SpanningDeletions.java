package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;

import java.util.List;


public class SpanningDeletions implements VariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        int deletions = 0;
        for (Integer offset : pileup.getOffsets() ) {
            if ( offset == -1 )
                deletions++;
        }
        return new Pair<String, String>("SpanningDeletions", String.format("%d", deletions));
    }

    public boolean useZeroQualityReads() { return false; }
}