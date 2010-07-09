package org.broadinstitute.sting.gatk.walkers.genotyper;

import static org.broadinstitute.sting.gatk.walkers.genotyper.BaseMismatchModel.*;

public class FourBaseProbabilitiesFactory {
    //private FourBaseProbabilitiesFactory() {} // cannot be instantiated

    public static BaseMismatchModel getBaseMismatchModel(final String name) {
        return valueOf(name);
    }

    public static BaseMismatchModel getBaseMismatchModel(final FourBaseProbabilities m) {
        if ( m instanceof OneStateErrorProbabilities)
            return ONE_STATE;
        else if ( m instanceof ThreeStateErrorProbabilities)
            return THREE_STATE;
        else if ( m instanceof EmpiricalSubstitutionProbabilities)
            return EMPIRICAL;

        throw new RuntimeException("Unexpected BaseMismatchModel " + m.getClass());

    }
    /**
     * General and correct way to create FourBaseLikelihood objects for arbitrary base mismatching models
     *
     * @param m  model
     * @param pl default platform
     * @return   new 4-base model
     */
    public static FourBaseProbabilities makeFourBaseLikelihoods(BaseMismatchModel m,
                                                              EmpiricalSubstitutionProbabilities.SequencerPlatform pl ) {
        switch ( m ) {
            case ONE_STATE: return new OneStateErrorProbabilities();
            case THREE_STATE: return new ThreeStateErrorProbabilities();
            case EMPIRICAL: return new EmpiricalSubstitutionProbabilities(pl);
            default: throw new RuntimeException("Unexpected BaseMismatchModel " + m);
        }
    }
}