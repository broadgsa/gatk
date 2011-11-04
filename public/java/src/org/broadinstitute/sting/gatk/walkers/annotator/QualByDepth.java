package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Variant confidence (given as (AB+BB)/AA from the PLs) / unfiltered depth.
 *
 * Low scores are indicative of false positive calls and artifacts.
 */
public class QualByDepth extends InfoFieldAnnotation implements StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        int depth = 0;

        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {

            // we care only about variant calls with likelihoods
            if ( genotype.getValue().isHomRef() )
                continue;

            AlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context == null )
                continue;

            depth += context.hasBasePileup() ? context.getBasePileup().depthOfCoverage() : context.getExtendedEventPileup().depthOfCoverage();
        }

        if ( depth == 0 )
            return null;

        double QD = 10.0 * vc.getNegLog10PError() / (double)depth;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", QD));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("QD"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Variant Confidence/Quality by Depth")); }

}