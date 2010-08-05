package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class GCContent implements InfoFieldAnnotation, ExperimentalAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        double content = computeGCContent(ref);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", content));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("GC"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("GC", 1, VCFHeaderLineType.Integer, "GC content within 20 bp +/- the variant")); }

    public boolean useZeroQualityReads() { return false; }

    private static double computeGCContent(ReferenceContext ref) {
        int gc = 0, at = 0;

        for ( byte base : ref.getBases() ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            if ( baseIndex == BaseUtils.gIndex || baseIndex == BaseUtils.cIndex )
                gc++;
            else if ( baseIndex == BaseUtils.aIndex || baseIndex == BaseUtils.tIndex )
                at++;
            else
                ; // ignore
        }

        int sum = gc + at;
        return (100.0*gc) / (sum == 0 ? 1 : sum);
     }
}