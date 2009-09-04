package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.refdata.RodGeliText;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;

public class VECAlleleBalance extends RatioFilter {

    private double lowThreshold = 0.25, highThreshold = 0.75;
    private double ratio;

    public VECAlleleBalance() {
        super("Allele Balance Ratio", VECAlleleBalance.class, Tail.TwoTailed);
    }

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            lowThreshold = Double.valueOf(argPieces[0]);
            highThreshold = Double.valueOf(argPieces[1]);
            if ( argPieces.length > 2 )
                minGenotypeConfidenceToTest = Double.valueOf(argPieces[2]);
        }
        setLowThreshold(lowThreshold);
        setHighThreshold(highThreshold);
    }

    /**
     * Return the count of bases matching the major (first) and minor (second) alleles as a pair.
     *
     */
    protected Pair<Integer, Integer> scoreVariant(char ref, ReadBackedPileup pileup, RodGeliText variant) {
        final String genotype = variant.getBestGenotype();
        final String bases = pileup.getBases();

        if ( genotype.length() > 2 )
            throw new IllegalArgumentException(String.format("Can only handle diploid genotypes: %s", genotype));


        char a = genotype.toUpperCase().charAt(0);
        char b = genotype.toUpperCase().charAt(1);
        int aCount = Utils.countOccurrences(a, bases.toUpperCase());
        int bCount = Utils.countOccurrences(b, bases.toUpperCase());

        int refCount = a == ref ? aCount : bCount;
        int altCount = a == ref ? bCount : aCount;

        ratio = (double)refCount / (double)altCount;
        return new Pair<Integer, Integer>(refCount, altCount);
    }

    protected boolean excludeHetsOnly() { return true; }

    public boolean useZeroQualityReads() { return false; }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "AlleleBalance("+lowThreshold+","+highThreshold+")\tRefRatio";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + ratio;
    }
}
