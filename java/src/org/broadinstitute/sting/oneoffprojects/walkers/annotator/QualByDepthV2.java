package org.broadinstitute.sting.oneoffprojects.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.AnnotationByDepth;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A Qual By Depth calculation adjusted to account for allele counts by (n var samples)/(n var alleles) -- so an entirely
 * homozygous variant receives a penalty of 1/2, while entirely het receives a (multiplicative) penalty of 1 (so no penalty)
 * This does not necessarily work well in the case of non-confident genotypes (could over or under penalize)
 */
public class QualByDepthV2 extends AnnotationByDepth implements ExperimentalAnnotation {
    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        //double QbyD = genotypeQualByDepth(genotypes, stratifiedContexts);
        int qDepth = annotationByVariantDepth(genotypes, stratifiedContexts);
        if ( qDepth == 0 )
            return null;

        double QbyD = hetHomAdjustment(vc) * 10.0 * vc.getNegLog10PError() / (double)qDepth;
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", QbyD));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("QD2"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth")); }

    public double hetHomAdjustment(VariantContext vc) {
        int variantSamples = 0;
        int variantAlleles = 0;
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            if ( ! g.isFiltered() ) {
                if ( g.isHet() ) {
                    variantSamples++;
                    variantAlleles++;
                } else if ( g.isHomVar() ) {
                    variantSamples++;
                    variantAlleles += 2;
                }
            }
        }

        return (0.0+variantSamples)/(0.0+variantAlleles);
    }
}
