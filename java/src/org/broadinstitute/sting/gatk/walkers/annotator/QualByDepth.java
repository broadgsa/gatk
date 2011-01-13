package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.Utils;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;


public class QualByDepth extends AnnotationByDepth implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        double qual = 0.0;
        int depth = 0;

        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {

            // we care only about variant calls with likelihoods
            if ( genotype.getValue().isHomRef() )
                continue;

            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context == null )
                continue;

            depth += context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();

            if ( genotype.getValue().hasLikelihoods() ) {
                GenotypeLikelihoods GLs = genotype.getValue().getLikelihoods();
                qual += 10.0 * getQual(GLs.getAsVector());
            }
        }

        if ( depth == 0 )
            return null;

        if ( qual == 0.0 )
            qual = 10.0 * vc.getNegLog10PError();

        double sumGLbyD = qual / (double)depth;

        int qDepth = annotationByVariantDepth(genotypes, stratifiedContexts);
        double QD = 10.0 * vc.getNegLog10PError() / (double)qDepth;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", QD));
        map.put(getKeyNames().get(1), String.format("%.2f", sumGLbyD));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("QD", "sumGLbyD"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth")); }

    private double getQual(double[] GLs) {
        if ( GLs == null )
            return 0.0;        

        // normalize so that we don't have precision issues
        double[] adjustedLikelihoods = new double[GLs.length];
        double maxValue = Utils.findMaxEntry(GLs);
        for (int i = 0; i < GLs.length; i++)
            adjustedLikelihoods[i] = GLs[i] - maxValue;

        // AB + BB (in real space)
        double variantWeight = Math.pow(10, adjustedLikelihoods[1]) + Math.pow(10, adjustedLikelihoods[2]);
        // (AB + BB) / AA (in log space)
        return Math.log10(variantWeight) - adjustedLikelihoods[0];
    }
 }