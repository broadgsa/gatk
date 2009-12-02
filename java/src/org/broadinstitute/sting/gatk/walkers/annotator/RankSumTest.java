package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.List;
import java.util.ArrayList;


public class RankSumTest implements VariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation, List<Genotype> genotypes) {

        if ( genotypes.size() == 0 )
            return null;

        // this test doesn't make sense for homs
        Genotype genotype = VariantAnnotator.getFirstHetVariant(genotypes);
        if ( genotype == null )
            return null;

        ArrayList<Integer> refQuals = new ArrayList<Integer>();
        ArrayList<Integer> altQuals = new ArrayList<Integer>();

        if ( genotype instanceof ReadBacked && ((ReadBacked)genotype).getPileup() != null )
            fillQualsFromGenotypes(ref.getBase(), variation.getAlternativeBaseForSNP(), genotypes, refQuals, altQuals);
        else
            fillQualsFromPileup(ref.getBase(), variation.getAlternativeBaseForSNP(), pileup, refQuals, altQuals);

        WilcoxonRankSum wilcoxon = new WilcoxonRankSum();
        for ( Integer qual : altQuals )
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET1);
        for ( Integer qual : refQuals )
            wilcoxon.addObservation((double)qual, WilcoxonRankSum.WILCOXON_SET.SET2);

        double pvalue = wilcoxon.getTwoTailedPValue();
        if ( MathUtils.compareDoubles(pvalue, 0.0) == 0 )
            return null;

        return String.format("%.1f", QualityUtils.phredScaleErrorRate(pvalue));
    }

    public String getKeyName() { return "RankSum"; }

    public String getDescription() { return "RankSum,1,Float,\"Phred-scaled p-value From Wilcoxon Rank Sum Test of Het Vs. Ref Base Qualities\""; }

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

    private void fillQualsFromGenotypes(char ref, char alt, List<Genotype> genotypes, List<Integer> refQuals, List<Integer> altQuals) {
        // accumulate quals
        for ( Genotype g : genotypes ) {
            if ( !(g instanceof ReadBacked) )
                continue;

            ReadBackedPileup pileup = ((ReadBacked)g).getPileup();
            if ( pileup == null )
                continue;

            fillQualsFromPileup(ref, alt, pileup, refQuals, altQuals);
        }
    }

    public boolean useZeroQualityReads() { return false; }
}