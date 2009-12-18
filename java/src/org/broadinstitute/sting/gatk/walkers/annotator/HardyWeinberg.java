package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.List;
import java.util.Map;


public class HardyWeinberg implements VariantAnnotation {

    private static final int MIN_SAMPLES = 10;
    private static final int MIN_GENOTYPE_QUALITY = 10;
    private static final int MIN_NEG_LOG10_PERROR = MIN_GENOTYPE_QUALITY / 10;

    public String annotate(ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {

        if ( !(variation instanceof VariantBackedByGenotype) )
            return null;
        final List<Genotype> genotypes = ((VariantBackedByGenotype)variation).getGenotypes();
        if ( genotypes == null || genotypes.size() < MIN_SAMPLES )
            return null;

        int refCount = 0;
        int hetCount = 0;
        int homCount = 0;
        for ( Genotype genotype : genotypes ) {
            if ( genotype.isNoCall() )
                continue;
            if ( genotype.getNegLog10PError() < MIN_NEG_LOG10_PERROR )
                continue;

            if ( !genotype.isVariant(ref.getBase()) )
                refCount++;
            else if ( genotype.isHet() )
                hetCount++;
            else
                homCount++;
        }

        double pvalue = HardyWeinbergCalculation.hwCalculate(refCount, hetCount, homCount);
        System.out.println(refCount + " " + hetCount + " " + homCount + " " + pvalue);
        return String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue));
    }

    public String getKeyName() { return "HW"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("HW", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Phred-scaled p-value for Hardy-Weinberg violation"); }
}