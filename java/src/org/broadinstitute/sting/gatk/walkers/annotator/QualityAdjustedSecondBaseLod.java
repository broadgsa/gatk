package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Dec 19, 2009
 * Time: 1:02:09 AM
 * To change this template use File | Settings | File Templates.
 */
public class QualityAdjustedSecondBaseLod implements VariantAnnotation {
    private final String KEY_NAME = "Qual_Adjusted_2blod";
    private final double CHI_LOD_MAX = -1000.0;
    private final SecondBaseSkew skewCalc = new SecondBaseSkew();
    private final double log10e = Math.log10(Math.E);
    private final double log10half = Math.log10(1.0/2);

    public String getKeyName() { return KEY_NAME; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(KEY_NAME, 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Adjusted residual quality based on second-base skew"); }
    
    public String annotate( RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> contexts, Variation variant) {
	String chi = skewCalc.annotate(tracker, ref,contexts,variant);
	if ( chi == null )
	    return null;
	double chi_square = Double.valueOf(chi);
        double chi_loglik = chi_square <= 0.0 ? 0.0 : Math.max(-(chi_square/2.0)*log10e + log10half,CHI_LOD_MAX); // cap it...
        return String.format("%f", 10*(variant.getNegLog10PError() + chi_loglik));
    }
}
