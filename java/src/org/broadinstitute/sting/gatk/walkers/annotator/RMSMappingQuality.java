package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.*;


public class RMSMappingQuality implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        ArrayList<Integer> qualities = new ArrayList<Integer>();
        for ( String sample : stratifiedContexts.keySet() ) {
            ReadBackedPileup pileup = stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for (PileupElement p : pileup )
                qualities.add(p.getRead().getMappingQuality());
        }
        int[] quals = new int[qualities.size()];
        int index = 0;
        for ( Integer i : qualities )
            quals[index++] = i;
        double rms = MathUtils.rms(quals);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", rms));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFRecord.RMS_MAPPING_QUALITY_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFInfoHeaderLine.INFO_TYPE.Float, "RMS Mapping Quality")); }
}