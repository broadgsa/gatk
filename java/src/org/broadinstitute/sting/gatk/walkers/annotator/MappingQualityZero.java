package org.broadinstitute.sting.gatk.walkers.annotator;

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


public class MappingQualityZero implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        int mq0 = 0;
        for ( String sample : stratifiedContexts.keySet() ) {
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for (PileupElement p : pileup ) {
                if ( p.getMappingQual() == 0 )
                    mq0++;
            }
        }
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%d", mq0));
        return map;
    }

    public String getKeyName() { return "MQ0"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Total Mapping Quality Zero Reads"); }
}