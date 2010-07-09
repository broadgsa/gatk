package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class MappingQualityZero implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int mq0 = 0;
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : stratifiedContexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for (PileupElement p : pileup ) {
                if ( p.getMappingQual() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", mq0));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("MQ0"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Integer, "Total Mapping Quality Zero Reads")); }
}