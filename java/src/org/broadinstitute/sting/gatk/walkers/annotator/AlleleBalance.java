package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.List;
import java.util.ArrayList;


public class AlleleBalance implements VariantAnnotation {

    public Pair<String, String> annotate(ReferenceContext ref, ReadBackedPileup pileup, List<Genotype> genotypes) {
        double ratio = -1;

        if ( genotypes.size() > 0 ) {
            Genotype g = genotypes.get(0);
            if ( g instanceof ReadBacked && g instanceof PosteriorsBacked ) {
                Pair<Double, Integer> weightedBalance = computeWeightedBalance(ref.getBase(), genotypes, pileup);
                if ( weightedBalance.second == 0 )
                    return null;
                ratio = weightedBalance.first;
            } else {
                // this test doesn't make sense for homs
                Genotype genotype = VariantAnnotator.getFirstHetVariant(genotypes);
                if ( genotype == null )
                    return null;

                final String genotypeStr = genotype.getBases().toUpperCase();
                if ( genotypeStr.length() != 2 )
                    return null;

                final String bases = pileup.getBasesAsString().toUpperCase();
                if ( bases.length() == 0 )
                    return null;

                ratio = computeSingleBalance(ref.getBase(), genotypeStr, bases);
            }
        }

        return new Pair<String, String>("AlleleBalance", String.format("%.2f", ratio));
    }

    private double computeSingleBalance(char ref, final String genotypeStr, final String bases) {

        char a = genotypeStr.charAt(0);
        char b = genotypeStr.charAt(1);
        int aCount = Utils.countOccurrences(a, bases);
        int bCount = Utils.countOccurrences(b, bases);

        int refCount = a == ref ? aCount : bCount;
        int altCount = a == ref ? bCount : aCount;

        double ratio = (double)refCount / (double)(refCount + altCount);
        return ratio;
    }

    // returns the ratio and then number of points which comprise it
    private Pair<Double, Integer> computeWeightedBalance(char ref, List<Genotype> genotypes, ReadBackedPileup pileup) {

        ArrayList<Double> refBalances = new ArrayList<Double>();
        ArrayList<Double> weights = new ArrayList<Double>();

        // accumulate ratios and weights
        for ( Genotype g : genotypes ) {

            if ( !(g instanceof ReadBacked) || !(g instanceof PosteriorsBacked) )
                continue;

            final String genotypeStr = g.getBases().toUpperCase();
            if ( genotypeStr.length() != 2 )
                continue;

            DiploidGenotype bestGenotype = DiploidGenotype.unorderedValueOf(genotypeStr);

            // we care only about het ref calls
            if ( !bestGenotype.isHetRef(ref) )
                continue;

            char a = genotypeStr.charAt(0);
            char b = genotypeStr.charAt(1);
            char altBase = a != ref ? a : b;

            // get the base counts at this pileup (minus deletions)
            ReadBackedPileup myPileup = ((ReadBacked)g).getPileup();

            // if the pileup is null, we'll just have to use the full pileup (it's better than nothing)
            if ( myPileup == null )
                myPileup = pileup;

            int[] counts = myPileup.getBasePileupAsCounts();
            int refCount = counts[BaseUtils.simpleBaseToBaseIndex(ref)];
            int altCount = counts[BaseUtils.simpleBaseToBaseIndex(altBase)];

            double[] posteriors = ((PosteriorsBacked)g).getPosteriors();
            posteriors = MathUtils.normalizeFromLog10(posteriors);
            weights.add(posteriors[bestGenotype.ordinal()]);
            refBalances.add((double)refCount / (double)(refCount + altCount));
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
                ratio += weights.get(i) * refBalances.get(i);
        }

        return new Pair<Double, Integer>(ratio, weights.size());
    }

    public boolean useZeroQualityReads() { return false; }
}
