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

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Dec 14, 2011
 *
 * Useful helper class to communicate the results of the allele frequency calculation
 */
public class AlleleFrequencyCalculationResult {

    // note that the cell at position zero in the likelihoods/posteriors array is actually probability of non-ref (since it's marginalized over all alleles)
    final double[][] log10AlleleFrequencyLikelihoods;
    final double[][] log10AlleleFrequencyPosteriors;

    double log10LikelihoodOfAFzero = 0.0;
    double log10PosteriorOfAFzero = 0.0;

    AlleleFrequencyCalculationResult(int maxAltAlleles, int numChr) {
        log10AlleleFrequencyLikelihoods = new double[maxAltAlleles][numChr+1];
        log10AlleleFrequencyPosteriors = new double[maxAltAlleles][numChr+1];
    }
}