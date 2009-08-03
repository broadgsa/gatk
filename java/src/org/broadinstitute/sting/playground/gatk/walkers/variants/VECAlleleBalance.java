package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;

public class VECAlleleBalance implements VariantExclusionCriterion {  //extends RatioFilter {
    //final private static GenotypeFeatureData.Tail tail = GenotypeFeatureData.Tail.TwoTailed;
    private boolean exclude;
    private double lowThreshold, highThreshold, ratio;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            lowThreshold = Double.valueOf(argPieces[0]);
            highThreshold = Double.valueOf(argPieces[1]);
        }
    }

    /**
     * We can only filter out het calls with the AlleleBalance filter
     *
     * @param variant
     * @return
     */
    protected boolean applyToVariant(rodVariants variant) {
        String genotype = variant.getBestGenotype();
        return genotype.charAt(0) != genotype.charAt(1);
    }

    /**
     * Return the count of bases matching the major (first) and minor (second) alleles as a pair.
     *
     * @param ref
     * @param pileup
     * @param variant
     * @return
     */
    protected Pair<Integer, Integer> scoreVariant(char ref, ReadBackedPileup pileup, rodVariants variant) {
        final String genotype = variant.getBestGenotype();
        final String bases = pileup.getBases();

        if ( genotype.length() > 2 )
            throw new IllegalArgumentException(String.format("Can only handle diploid genotypes: %s", genotype));


        char a = genotype.toUpperCase().charAt(0);
        char b = genotype.toUpperCase().charAt(1);
        int aCount = Utils.countOccurrences(a, bases.toUpperCase());
        int bCount = Utils.countOccurrences(b, bases.toUpperCase());

        int major = Math.max(aCount, bCount);
        int minor = Math.min(aCount, bCount);

        return new Pair(major, minor);
    }

    public void compute(char ref, LocusContext context, rodVariants variant) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        Pair<Integer, Integer> counts = scoreVariant(ref, pileup, variant);
        ratio = (double)counts.first / (double)counts.second;
        exclude = ratio < lowThreshold || ratio > highThreshold;
    }

    public boolean useZeroQualityReads() { return false; }

    public boolean isExcludable() {
        return exclude;
    }

    public String getStudyHeader() {
        return "AlleleBalance("+lowThreshold+","+highThreshold+")\tMajorMinorRatio";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + ratio;
    }
}
