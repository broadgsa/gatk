package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;


public abstract class RankSumTest implements InfoFieldAnnotation, WorkInProgressAnnotation {
    private final static boolean DEBUG = false;
    private static final double minPValue = 1e-10;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {

        if ( !vc.isBiallelic() || !vc.isSNP() )
            return null;
        
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        ArrayList<Integer> refQuals = new ArrayList<Integer>();
        ArrayList<Integer> altQuals = new ArrayList<Integer>();

        for ( Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
            // we care only about het calls
            if ( genotype.getValue().isHet() ) {
                StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
                if ( context == null )
                    continue;

                fillQualsFromPileup(ref.getBaseAsChar(), vc.getAlternateAllele(0).toString().charAt(0), context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup(), refQuals, altQuals);
            }
        }

        WilcoxonRankSum wilcoxon = new WilcoxonRankSum();
        for ( Integer qual : altQuals )
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET1);
        for ( Integer qual : refQuals )
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET2);

        // for R debugging
        if ( DEBUG ) {
            wilcoxon.DEBUG = DEBUG;
            System.out.printf("%s%n", ref.getLocus());
            System.out.printf("alt <- c(%s)%n", Utils.join(",", altQuals));
            System.out.printf("ref <- c(%s)%n", Utils.join(",", refQuals));
        }

        // we are testing these set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        double pvalue = wilcoxon.getPValue(WilcoxonRankSum.WILCOXON_H0.SMALLER_SET_LT);
        if ( MathUtils.compareDoubles(pvalue, -1.0) == 0 )
            return null;

        // deal with precision issues
        if ( pvalue < minPValue )
            pvalue = minPValue;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyName(), String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue)));
        return map;
    }

    protected abstract void fillQualsFromPileup(char ref, char alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals);
}