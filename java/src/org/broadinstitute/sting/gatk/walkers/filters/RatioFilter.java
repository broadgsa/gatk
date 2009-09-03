package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.utils.*;
import org.apache.log4j.Logger;
import cern.jet.math.Arithmetic;


public abstract class RatioFilter implements VariantExclusionCriterion {
    private static final double defaultMinGenotypeConfidenceToTest = 5.0;  // TODO -- must be replaced with true confidence scoore, right now assumes LOD
    private static final int minDepthOfCoverage = 25;                   // TODO -- must be replaced with a proper probability calculation
    protected double minGenotypeConfidenceToTest = defaultMinGenotypeConfidenceToTest;

    protected double pvalueLimit = -1;
    protected Logger logger = null;    // Logger.getLogger(RatioFilter.class);
    protected String name = null;

    protected enum Tail { LeftTailed, RightTailed, TwoTailed }
    protected Tail tail = null;
    protected double lowThreshold = -1;
    protected double highThreshold = -1;

    protected boolean exclude = false;

    public RatioFilter(final String name, Class myClass, Tail tail ) {
        this.name = name;
        this.tail = tail;
        logger = Logger.getLogger(myClass);
    }

    protected void setLowThreshold(double threshold) {
        lowThreshold = threshold;
    }

    protected void setHighThreshold(double threshold) {
        highThreshold = threshold;
    }

    protected abstract Pair<Integer, Integer> scoreVariant(char ref, ReadBackedPileup pileup, rodVariants variant);
    protected abstract boolean excludeHetsOnly();

    public boolean useZeroQualityReads() { return false; }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        rodVariants variant = context.getVariant();
        char ref = context.getReferenceContext().getBase();

        ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getAlignmentContext(useZeroQualityReads()));
        Pair<Integer, Integer> counts = scoreVariant(ref, pileup, variant);

        boolean highGenotypeConfidence = variant.getConsensusConfidence() > minGenotypeConfidenceToTest;
        boolean excludable = !excludeHetsOnly() || GenotypeUtils.isHet(variant);
        exclude = excludable && highGenotypeConfidence && pointEstimateExclude(counts);
        //
        // for printing only
        //
        int n = counts.first + counts.second;
        double value = counts.first / (1.0 * counts.first + counts.second);
        logger.info(String.format("%s: counts1=%d (%.2f), counts2=%d (%.2f), n=%d, value=%f, exclude=%b, location=%s, bases=%s",
                name, counts.first, counts.first / (0.01 * n), counts.second, counts.second / (0.01 * n), n,
                value, exclude, variant.getLocation(), pileup.getBases()));
    }

    private final static double SEARCH_INCREMENT = 0.01;
    private final static double integralPValueThreshold = 0.05;

    // TODO - this whole calculation needs to be redone correctly
    private boolean pointEstimateExclude(Pair<Integer, Integer> counts) {
        if ( counts.first + counts.second < minDepthOfCoverage )
            return false;

        int n = counts.first + counts.second;
        double ratio = counts.first.doubleValue() / (double)n;
        return !passesThreshold(ratio);
    }

    private boolean integralExclude(Pair<Integer, Integer> counts) {
        double sumExclude = 0.0, sumP = 0.0;
        int n = counts.first + counts.second;
        for ( double r = 0.0; r <= 1.0; r += SEARCH_INCREMENT ) {
            double p = MathUtils.binomialProbability(counts.first, n, r);
            sumP += p;
            boolean exclude = ! passesThreshold(r);
            sumExclude += p * (exclude ? 1.0 : 0.0);

            //System.out.printf("integral: k=%d, n=%d, r=%f, p=%f, sumP = %f, exclude=%b | sum=%f, percentExcluded=%f%n",
            //        counts.first, n, r, p, sumP, exclude, sumExclude, sumExclude / sumP);
        }

        double percentExcluded = sumExclude / sumP;
        return 1 - percentExcluded <= integralPValueThreshold ;
    }

    private boolean passesThreshold(double value) {
        switch ( tail ) {
            case LeftTailed:
                return value >= lowThreshold;
            case RightTailed:
                return value <= highThreshold;
            case TwoTailed:
                return value >= lowThreshold && value < highThreshold;
            default:
                return true;
        }
    }
}