package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.ReadBacked;
import org.broadinstitute.sting.utils.genotype.PosteriorsBacked;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import java.util.List;
import java.util.ArrayList;

public class OnOffGenotype implements VariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {

        if ( genotypes.size() == 0 )
            return null;

        double ratio;
        Genotype g = genotypes.get(0);
        if ( g instanceof ReadBacked && g instanceof PosteriorsBacked) {
            Pair<Double, Integer> weightedBalance = computeWeightedBalance(ref.getBase(), genotypes, pileup);
            if ( weightedBalance.second == 0 )
                return null;
            ratio = weightedBalance.first;
        } else {
            Genotype genotype = VariantAnnotator.getFirstVariant(ref.getBase(), genotypes);
            if ( genotype == null )
                return null;

            final String genotypeStr = genotype.getBases().toUpperCase();
            if ( genotypeStr.length() != 2 )
                return null;

            final String bases = pileup.getBasesAsString().toUpperCase();
            if ( bases.length() == 0 )
                return null;

            ratio = computeSingleBalance(genotypeStr, bases);
        }

        return new Pair<String, String>("OnOffGenotype", String.format("%.2f", ratio));
    }

    private double computeSingleBalance(final String genotypeStr, final String bases) {

        int on = 0, off = 0;
        for ( char base : BaseUtils.BASES ) {
            int count = BasicPileup.countBase(base, bases);
            if ( Utils.countOccurrences(base, genotypeStr) > 0 )
                on += count;
            else
                off += count;
        }

        double ratio = (double)on / (double)(on + off);
        return ratio;
    }

    private Pair<Double, Integer> computeWeightedBalance(char ref, List<Genotype> genotypes, ReadBackedPileup pileup) {

        ArrayList<Double> onOffBalances = new ArrayList<Double>();
        ArrayList<Double> weights = new ArrayList<Double>();

        // accumulate ratios and weights
        for ( Genotype g : genotypes ) {

            if ( !(g instanceof ReadBacked) || !(g instanceof PosteriorsBacked) )
                continue;

            final String genotypeStr = g.getBases().toUpperCase();
            if ( genotypeStr.length() != 2 )
                continue;

            DiploidGenotype bestGenotype = DiploidGenotype.unorderedValueOf(genotypeStr);

            // we care only about non-ref calls
            if ( bestGenotype.isHomRef(ref) )
                continue;

            char a = genotypeStr.charAt(0);
            char b = genotypeStr.charAt(1);

            // get the base counts at this pileup (minus deletions)
            ReadBackedPileup myPileup = ((ReadBacked)g).getPileup();

            // if the pileup is null, we'll just have to use the full pileup (it's better than nothing)
            if ( myPileup == null )
                myPileup = pileup;

            int[] counts = myPileup.getBasePileupAsCounts();
            int onCount = counts[BaseUtils.simpleBaseToBaseIndex(a)];
            if ( a != b )
                onCount += counts[BaseUtils.simpleBaseToBaseIndex(b)];
            int totalCount = 0;
            for (int i = 0; i < counts.length; i++)
                totalCount += counts[i];

            // sanity check
            if ( totalCount == 0 )
                continue;

            double[] posteriors = ((PosteriorsBacked)g).getPosteriors();
            posteriors = MathUtils.normalizeFromLog10(posteriors);
            double weight = posteriors[bestGenotype.ordinal()];

            // sanity check
            if ( MathUtils.compareDoubles(weight, 0.0) == 0 )
                continue;

            weights.add(weight);
            onOffBalances.add((double)onCount / (double)totalCount);
        }

        double ratio = 0.0;

        if ( weights.size() > 0 ) {
            // normalize the weights
            double sum = 0.0;
            for (int i = 0; i < weights.size(); i++)
                sum += weights.get(i);
            for (int i = 0; i < weights.size(); i++)
                weights.set(i, weights.get(i) / sum);

            // calculate total weighted ratios
            for (int i = 0; i < weights.size(); i++)
                ratio += weights.get(i) * onOffBalances.get(i);
        }

        return new Pair<Double, Integer>(ratio, weights.size());
    }

    public boolean useZeroQualityReads() { return false; }
}
