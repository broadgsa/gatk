package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;


public class SpanningDeletions extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        int deletions = 0;
        int depth = 0;
        for ( String sample : stratifiedContexts.keySet() ) {
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.MQ0FREE).getPileup();
            deletions += pileup.getNumberOfDeletions();
            depth += pileup.size();
        }
        return String.format("%.2f", (double)deletions/(double)depth);
    }

    public String getKeyName() { return "Dels"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("Dels", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Fraction of Reads Containing Spanning Deletions"); }
}