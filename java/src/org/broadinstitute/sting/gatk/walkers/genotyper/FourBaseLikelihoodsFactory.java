/*
 * Copyright (c) 2010.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.genotyper;

import static org.broadinstitute.sting.gatk.walkers.genotyper.BaseMismatchModel.*;

public class FourBaseLikelihoodsFactory {
    //private FourBaseProbabilitiesFactory() {} // cannot be instantiated

    public static BaseMismatchModel getBaseMismatchModel(final String name) {
        return valueOf(name);
    }

    public static BaseMismatchModel getBaseMismatchModel(final FourBaseLikelihoods m) {
        if ( m instanceof ThreeStateErrorProbabilities)
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
    public static FourBaseLikelihoods makeFourBaseLikelihoods(BaseMismatchModel m,
                                                              EmpiricalSubstitutionProbabilities.SequencerPlatform pl ) {
        switch ( m ) {
            case THREE_STATE: return new ThreeStateErrorProbabilities();
            case EMPIRICAL: return new EmpiricalSubstitutionProbabilities(pl);
            default: throw new RuntimeException("Unexpected BaseMismatchModel " + m);
        }
    }
}