package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;


public class QualByDepth implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        //double QbyD = genotypeQualByDepth(genotypes, stratifiedContexts);
        int qDepth = variationQualByDepth(genotypes, stratifiedContexts);
        if ( qDepth == 0 )
            return null;

        double QbyD = 10.0 * vc.getNegLog10PError() / (double)qDepth;
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%.2f", QbyD));
        return map;
    }

    public String getKeyName() { return "QD"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Variant Confidence/Quality by Depth"); }

    private int variationQualByDepth(final Map<String, Genotype> genotypes, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
        int depth = 0;
        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {

            // we care only about variant calls
            if ( genotype.getValue().isHomRef() )
                continue;

            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context != null )
                depth += context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();
        }

        return depth;
    }

//    private double genotypeQualByDepth(final Map<String, Genotype> genotypes, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
//        ArrayList<Double> qualsByDepth = new ArrayList<Double>();
//        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
//
//            // we care only about variant calls
//            if ( genotype.getValue().isHomRef() )
//                continue;
//
//            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
//            if ( context == null )
//                continue;
//
//            qualsByDepth.add(genotype.getValue().getNegLog10PError() / context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size());
//        }
//
//        double mean = 0.0;
//        for ( Double value : qualsByDepth )
//            mean += value;
//        mean /= qualsByDepth.size();
//
//        return mean;
//    }
}