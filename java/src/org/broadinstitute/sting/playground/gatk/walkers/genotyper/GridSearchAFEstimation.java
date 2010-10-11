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

package org.broadinstitute.sting.playground.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.*;
import java.io.PrintStream;

public class GridSearchAFEstimation extends AlleleFrequencyCalculationModel {

    // for use in optimizing the P(D|AF) calculations:
    // how much off from the max likelihoods do we need to be before we can quit calculating?
    protected static final double LOG10_OPTIMIZATION_EPSILON = 8.0;

    protected GridSearchAFEstimation(int N, Logger logger, PrintStream verboseWriter) {
        super(N, logger, verboseWriter);
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, BiallelicGenotypeLikelihoods> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors) {

        initializeAFMatrix(GLs);

        // first, calculate for AF=0 (no change to matrix)
        log10AlleleFrequencyPosteriors[0] = AFMatrix.getLikelihoodsOfFrequency() + log10AlleleFrequencyPriors[0];
        double maxLikelihoodSeen = log10AlleleFrequencyPosteriors[0];

        int maxAlleleFrequencyToTest = AFMatrix.getSamples().size() * 2;

        // for each minor allele frequency, calculate log10PofDgivenAFi
        for (int i = 1; i <= maxAlleleFrequencyToTest; i++) {
            // add one more alternate allele
            AFMatrix.incrementFrequency();

            // calculate new likelihoods
            log10AlleleFrequencyPosteriors[i] = AFMatrix.getLikelihoodsOfFrequency() + log10AlleleFrequencyPriors[i];

            // an optimization to speed up the calculation: if we are beyond the local maximum such
            //  that subsequent likelihoods won't factor into the confidence score, just quit
            if ( maxLikelihoodSeen - log10AlleleFrequencyPosteriors[i] > LOG10_OPTIMIZATION_EPSILON )
                return;

            if ( log10AlleleFrequencyPosteriors[i] > maxLikelihoodSeen )
                maxLikelihoodSeen = log10AlleleFrequencyPosteriors[i];
        }
    }

    /**
     * Overrides the super class
     * @param contexts             alignment contexts
     * @param GLs                  genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     *
     * @return calls
     */
    public Map<String, Genotype> assignGenotypes(Map<String, StratifiedAlignmentContext> contexts,
                                                 Map<String, BiallelicGenotypeLikelihoods> GLs,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int AFofMaxLikelihood) {
        return generateCalls(contexts, GLs, AFofMaxLikelihood);
    }
}