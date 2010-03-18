package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;
import java.util.HashMap;

public class QualityAdjustedSecondBaseLod implements InfoFieldAnnotation, ExperimentalAnnotation {
    private final String KEY_NAME = "Qual_Adjusted_2blod";
    private final double CHI_LOD_MAX = -1000.0;
    private final SecondBaseSkew skewCalc = new SecondBaseSkew();
    private final double log10e = Math.log10(Math.E);
    private final double log10half = Math.log10(1.0/2);

    public String getKeyName() { return KEY_NAME; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(KEY_NAME, 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Adjusted residual quality based on second-base skew"); }
    
    public Map<String, Object> annotate( RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> contexts, VariantContext vc) {
        String chi = skewCalc.getAnnotation(ref, contexts, vc);
        if ( chi == null )
            return null;
        double chi_square = Double.valueOf(chi);
        double chi_loglik = chi_square <= 0.0 ? 0.0 : Math.max(-(chi_square/2.0)*log10e + log10half,CHI_LOD_MAX); // cap it...
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%f", 10*(vc.getNegLog10PError() + chi_loglik)));
        return map;
    }
}
