package org.broadinstitute.sting.gatk.walkers.genotyper;

import static org.broadinstitute.sting.gatk.walkers.genotyper.BaseMismatchModel.*;

public class GenotypeLikelihoodsFactory {
    //private GenotypeLikelihoodsFactory() {} // cannot be instantiated

    public static BaseMismatchModel getBaseMismatchModel(final String name) {
        BaseMismatchModel m = valueOf(name);
        if ( m == null )
            throw new RuntimeException("Unexpected BaseMismatchModel " + name);
        else
            return m;
    }

    /**
     * General and correct way to create GenotypeLikelihood objects for arbitrary base mismatching models
     *
     * @param m
     * @param priors
     * @param pl
     * @return
     */
    public static GenotypeLikelihoods makeGenotypeLikelihoods(BaseMismatchModel m,
                                                              DiploidGenotypePriors priors,
                                                              EmpiricalSubstitutionGenotypeLikelihoods.SequencerPlatform pl ) {
        switch ( m ) {
            case ONE_STATE: return new OneStateErrorGenotypeLikelihoods(priors);
            case THREE_STATE: return new ThreeStateErrorGenotypeLikelihoods(priors);
            case EMPIRICAL: return new EmpiricalSubstitutionGenotypeLikelihoods(priors, pl);
            default: throw new RuntimeException("Unexpected BaseMismatchModel " + m);
        }
    }
}