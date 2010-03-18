package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;
import java.util.HashMap;


public class AlleleBalance implements InfoFieldAnnotation, StandardAnnotation {

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( !vc.isBiallelic() || !vc.isSNP() )
            return null;
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        double ratio = 0.0;
        double totalWeights = 0.0;
        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
            // we care only about het calls
            if ( !genotype.getValue().isHet() )
                continue;

            StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context == null )
                continue;

            final String bases = new String(context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().getBases()).toUpperCase();
            if ( bases.length() == 0 )
                return null;

            char refChr = vc.getReference().toString().charAt(0);
            char altChr = vc.getAlternateAllele(0).toString().charAt(0);

            int refCount = Utils.countOccurrences(refChr, bases);
            int altCount = Utils.countOccurrences(altChr, bases);

            // sanity check
            if ( refCount + altCount == 0 )
                continue;

            // weight the allele balance by genotype quality so that e.g. mis-called homs don't affect the ratio too much
            ratio += genotype.getValue().getNegLog10PError() * ((double)refCount / (double)(refCount + altCount));
            totalWeights += genotype.getValue().getNegLog10PError();
        }

        // make sure we had a het genotype
        if ( MathUtils.compareDoubles(totalWeights, 0.0) == 0 )
            return null;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%.2f", (ratio / totalWeights)));
        return map;
    }

    public String getKeyName() { return "AB"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("AB", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Balance for hets (ref/(ref+alt))"); }
}
