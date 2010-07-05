package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.genotype.HardyWeinbergCalculation;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class HardyWeinberg implements InfoFieldAnnotation, WorkInProgressAnnotation {

    private static final int MIN_SAMPLES = 10;
    private static final int MIN_GENOTYPE_QUALITY = 10;
    private static final int MIN_NEG_LOG10_PERROR = MIN_GENOTYPE_QUALITY / 10;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() < MIN_SAMPLES )
            return null;

        int refCount = 0;
        int hetCount = 0;
        int homCount = 0;
        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
            Genotype g = genotype.getValue();

            if ( g.isNoCall() )
                continue;

            // TODO - fix me:
            // Right now we just ignore genotypes that are not confident, but this throws off
            //  our HW ratios.  More analysis is needed to determine the right thing to do when
            //  the genotyper cannot decide whether a given sample is het or hom var.
            if ( g.getNegLog10PError() < MIN_NEG_LOG10_PERROR )
                continue;

            if ( g.isHomRef() )
                refCount++;
            else if ( g.isHet() )
                hetCount++;
            else
                homCount++;
        }

        if ( refCount + hetCount + homCount == 0)
            return null;

        double pvalue = HardyWeinbergCalculation.hwCalculate(refCount, hetCount, homCount);
        //System.out.println(refCount + " " + hetCount + " " + homCount + " " + pvalue);
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue)));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("HW"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("HW", 1, VCFHeaderLineType.Float, "Phred-scaled p-value for Hardy-Weinberg violation")); }
}