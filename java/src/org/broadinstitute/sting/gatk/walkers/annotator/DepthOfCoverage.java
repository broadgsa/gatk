package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class DepthOfCoverage implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        int depth = 0;
        for ( String sample : stratifiedContexts.keySet() )
            depth += stratifiedContexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%d", depth));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(VCFRecord.DEPTH_KEY); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Total Depth")); }
}
