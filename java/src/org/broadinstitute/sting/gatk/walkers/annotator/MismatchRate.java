package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.*;

import java.util.Map;


public class MismatchRate implements VariantAnnotation {

    public String annotate(ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        int mismatches = 0;
        int totalBases = 0;
        for ( String sample : stratifiedContexts.keySet() ) {
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.MQ0FREE).getPileup();
            Pair<Integer, Integer> counts = AlignmentUtils.mismatchesInRefWindow(pileup, ref, true);
            mismatches += counts.first;
            totalBases += counts.second;
        }

        // sanity check
        if ( totalBases == 0 )
            return null;

        return String.format("%.2f", (double)mismatches/(double)totalBases);
    }

    public String getKeyName() { return "MR"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("MR", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Mismatch Rate of Reads Spanning This Position"); }
}