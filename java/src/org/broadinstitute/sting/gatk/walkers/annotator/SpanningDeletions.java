package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class SpanningDeletions implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int deletions = 0;
        int depth = 0;
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : stratifiedContexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            deletions += pileup.getNumberOfDeletions();
            depth += pileup.size();
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", depth == 0 ? 0.0 : (double)deletions/(double)depth));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("Dels"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("Dels", 1, VCFHeaderLineType.Float, "Fraction of Reads Containing Spanning Deletions")); }
}