package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenotypeUtils;

public class VECAlleleBalance implements VariantExclusionCriterion {  //extends RatioFilter {
    private static final double defaultMinGenotypeConfidenceToTest = 5.0;  // TODO -- must be replaced with true confidence scoore, right now assumes LOD

    private boolean exclude;
    private double lowThreshold = 0.1, highThreshold = 0.85;
    private double minGenotypeConfidenceToTest = defaultMinGenotypeConfidenceToTest;
    private double ratio;

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            lowThreshold = Double.valueOf(argPieces[0]);
            highThreshold = Double.valueOf(argPieces[1]);
            if ( argPieces.length > 2 )
                minGenotypeConfidenceToTest = Double.valueOf(argPieces[2]);
        }
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

        int refCount = a == ref ? aCount : bCount;
        int altCount = a == ref ? bCount : aCount;

        return new Pair(refCount, altCount);
    }

    public void compute(char ref, AlignmentContext context, rodVariants variant) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        Pair<Integer, Integer> counts = scoreVariant(ref, pileup, variant);
        int n = counts.first + counts.second;
        ratio = counts.first.doubleValue() / (double)n;

        boolean highGenotypeConfidence = variant.getConsensusConfidence() > minGenotypeConfidenceToTest;
        boolean failsHetExpectation = GenotypeUtils.isHet(variant) && (ratio < lowThreshold || ratio > highThreshold);
        exclude = failsHetExpectation && highGenotypeConfidence;

//        if ( failsHetExpectation ) {
//            String header = highGenotypeConfidence ? "FILTER-HIGH-CONFIDENCE" : "PASS-LOW-CONFIDENCE";
//            System.out.printf("[%s] %s getConsensusConfidence() = %f > minGenotypeConfidenceToTest = %f exclude=%b  %s%n",
//                    header, variant.getBestGenotype(), variant.getConsensusConfidence(), minGenotypeConfidenceToTest, exclude,
//                    ! highGenotypeConfidence ? variant : "" );
//        }
    }

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
