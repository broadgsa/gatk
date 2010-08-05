package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;



public class SBByDepth extends AnnotationByDepth {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        if (!vc.hasAttribute(VCFConstants.STRAND_BIAS_KEY))
            return null;

        double sBias = Double.valueOf(vc.getAttributeAsString(VCFConstants.STRAND_BIAS_KEY));

        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        int sDepth = annotationByVariantDepth(genotypes, stratifiedContexts);
        if ( sDepth == 0 )
            return null;



        double SbyD = (-sBias / (double)sDepth);
        if (SbyD > 0)
            SbyD = Math.log10(SbyD);
        else
            SbyD = -1000;
        
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.2f", SbyD));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList("SBD"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Strand Bias by Depth")); }



}