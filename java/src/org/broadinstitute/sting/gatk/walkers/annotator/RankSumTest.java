package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;


public class RankSumTest implements VariantAnnotation {
    private final static boolean DEBUG = false;
    private static final double minPValue = 1e-10;

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {

        if ( !(variation instanceof VariantBackedByGenotype) )
            return null;
        final List<Genotype> genotypes = ((VariantBackedByGenotype)variation).getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        ArrayList<Integer> refQuals = new ArrayList<Integer>();
        ArrayList<Integer> altQuals = new ArrayList<Integer>();

        for ( Genotype genotype : genotypes ) {
            // we care only about het calls
            if ( genotype.isHet() ) {
                String sample = ((SampleBacked)genotype).getSampleName();
                StratifiedAlignmentContext context = stratifiedContexts.get(sample);
                if ( context == null )
                    continue;

                fillQualsFromPileup(ref.getBase(), variation.getAlternativeBaseForSNP(), context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup(), refQuals, altQuals);
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

        return String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue));
    }

    public String getKeyName() { return "RankSum"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("RankSum", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Phred-scaled p-value From Wilcoxon Rank Sum Test of Het Vs. Ref Base Qualities"); }

    private void fillQualsFromPileup(char ref, char alt, ReadBackedPileup pileup, List<Integer> refQuals, List<Integer> altQuals) {
        for ( PileupElement p : pileup ) {
            // ignore deletions
            if ( p.isDeletion() )
                continue;

            char base = (char)p.getBase();
            if ( base == ref )
                refQuals.add((int)p.getQual());
            else if ( base == alt )
                altQuals.add((int)p.getQual());
        }
    }
}