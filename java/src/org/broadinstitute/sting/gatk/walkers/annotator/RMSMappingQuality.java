package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.util.Map;
import java.util.ArrayList;


public class RMSMappingQuality extends StandardVariantAnnotation {

    public String annotate(ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        ArrayList<Integer> qualities = new ArrayList<Integer>();
        for ( String sample : stratifiedContexts.keySet() ) {
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getPileup();
            for (PileupElement p : pileup )
                qualities.add(p.getRead().getMappingQuality());
        }
        int[] quals = new int[qualities.size()];
        int index = 0;
        for ( Integer i : qualities )
            quals[index++] = i;
        double rms = MathUtils.rms(quals);
        return String.format("%.2f", rms);
    }

    public String getKeyName() { return VCFRecord.RMS_MAPPING_QUALITY_KEY; }

    public String getDescription() { return getKeyName() + ",1,Float,\"RMS Mapping Quality\""; }
}