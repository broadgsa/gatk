package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;


public abstract class RankSumTest implements InfoFieldAnnotation, ExperimentalAnnotation {
    private static final double minPValue = 1e-20;

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;
         
        if ( !vc.isBiallelic() || !vc.isSNP() )
            return null;
        
        final Map<String, Genotype> genotypes = vc.getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        final ArrayList<Integer> refQuals = new ArrayList<Integer>();
        final ArrayList<Integer> altQuals = new ArrayList<Integer>();

        for ( final Map.Entry<String, Genotype> genotype : genotypes.entrySet() ) {
            final StratifiedAlignmentContext context = stratifiedContexts.get(genotype.getKey());
            if ( context == null ) {
                continue;
            }
            fillQualsFromPileup(ref.getBase(), vc.getAlternateAllele(0).getBases()[0], context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup(), refQuals, altQuals);
        }

        final MannWhitneyU mannWhitneyU = new MannWhitneyU();
        for ( final Integer qual : altQuals ) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET1);
        }
        for ( final Integer qual : refQuals ) {
            mannWhitneyU.add(qual, MannWhitneyU.USet.SET2);
        }

        // we are testing that set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        final Pair<Double,Double> testResults = mannWhitneyU.runOneSidedTest( MannWhitneyU.USet.SET1 );
        double pvalue = testResults.second;

        /* // WilcoxonRankSum.java is slated for removal
        final WilcoxonRankSum wilcoxon = new WilcoxonRankSum();
        for ( final Integer qual : altQuals ) {
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET1);
        }
        for ( final Integer qual : refQuals ) {
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET2);
        }

        // we are testing these set1 (the alt bases) have lower quality scores than set2 (the ref bases)
        double pvalue = wilcoxon.getPValue(WilcoxonRankSum.WILCOXON_H0.SET1_LT_SET2);

        if ( MathUtils.compareDoubles(pvalue, -1.0) == 0 ) {
            pvalue = 1.0;
        }
        */

        // deal with precision issues
        if ( pvalue < minPValue ) {
            pvalue = minPValue;
        }

        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.3f", Math.abs(QualityUtils.phredScaleErrorRate(pvalue))));
        return map;
    }

    protected abstract void fillQualsFromPileup(byte ref, byte alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals);

    protected static boolean isUsableBase( final PileupElement p ) {
        return !( p.isDeletion() || p.getMappingQual() == 0 || ((int)p.getQual()) < 6 ); // need the unBAQed quality score here
    }
}