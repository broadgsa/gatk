package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.refdata.RodGeliText;
import org.broadinstitute.sting.utils.*;

import java.util.HashMap;

public class VECOnOffGenotypeRatio extends RatioFilter {
    private double threshold = 0.8;
    private double ratio;

    public VECOnOffGenotypeRatio() {
        super("On/Off Genotype Ratio", VECOnOffGenotypeRatio.class, Tail.LeftTailed);
    }

    public void initialize(HashMap<String,String> args) {
        if ( args.get("threshold") != null )
            threshold = Double.valueOf(args.get("threshold"));
        if ( args.get("confidence") != null )
            minGenotypeConfidenceToTest = Double.valueOf(args.get("confidence"));
        setLowThreshold(threshold);
    }

    /**
     * Return the counts of bases that are on (matching the bestGenotype) and off (not matching the
     * best genotype).  On are in the first field, off in the second.
     *
    */
    protected Pair<Integer, Integer> scoreVariant(char ref, ReadBackedPileup pileup, RodGeliText variant) {
        final String genotype = variant.getBestGenotype().toUpperCase();
        final String bases = pileup.getBases();

        if ( genotype.length() > 2 )
            throw new IllegalArgumentException(String.format("Can only handle diploid genotypes: %s", genotype));

        int on = 0, off = 0;

        for ( char base : BaseUtils.BASES ) {
            int count = BasicPileup.countBase(base, bases);
            if ( Utils.countOccurrences(base, genotype) > 0 )
                on += count;
            else
                off += count;
            //System.out.printf("count = %d, on=%d, off=%d for %c in %s%n", count, on, off, base, genotype);            
        }

        ratio = (double)on / (double)off;
        return new Pair<Integer, Integer>(on, off);
    }

    protected boolean excludeHetsOnly() { return false; }

    public double inclusionProbability() {
        return exclude ? 0.0 : 1.0;
    }

    public boolean useZeroQualityReads() { return false; }

    public String getStudyHeader() {
        return "OnOffGenotype("+threshold+")\tOnRatio";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + ratio;
    }
}