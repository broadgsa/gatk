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

import org.broadinstitute.sting.utils.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Dec 14, 2011
 *
 * Useful helper class to communicate the results of the allele frequency calculation
 */
public class AlleleFrequencyCalculationResult {

    // These variables are intended to contain the MLE and MAP (and their corresponding allele counts) of the site over all alternate alleles
    private double log10MLE;
    private double log10MAP;
    final private int[] alleleCountsOfMLE;
    final private int[] alleleCountsOfMAP;

    // The posteriors seen, not including that of AF=0
    // TODO -- better implementation needed here (see below)
    private ArrayList<Double> log10PosteriorMatrixValues = new ArrayList<Double>(100000);
    private Double log10PosteriorMatrixSum = null;

    // These variables are intended to contain the likelihood/posterior probability for the site's being monomorphic (i.e. AF=0 for all alternate alleles)
    private double log10LikelihoodOfAFzero;
    private double log10PosteriorOfAFzero;


    public AlleleFrequencyCalculationResult(final int maxAltAlleles) {
        alleleCountsOfMLE = new int[maxAltAlleles];
        alleleCountsOfMAP = new int[maxAltAlleles];
        reset();
    }

    public double getLog10MLE() {
        return log10MLE;
    }

    public double getLog10MAP() {
        return log10MAP;
    }

    public double getLog10PosteriorMatrixSum() {
        if ( log10PosteriorMatrixSum == null ) {
            // TODO -- we absolutely need a better implementation here as we don't want to store all values from the matrix in memory;
            // TODO --   will discuss with the team what the best option is
            final double[] tmp = new double[log10PosteriorMatrixValues.size()];
            for ( int i = 0; i < tmp.length; i++ )
                tmp[i] = log10PosteriorMatrixValues.get(i);
            log10PosteriorMatrixSum = MathUtils.log10sumLog10(tmp);
        }
        return log10PosteriorMatrixSum;
    }

    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE;
    }

    public int[] getAlleleCountsOfMAP() {
        return alleleCountsOfMAP;
    }

    public double getLog10LikelihoodOfAFzero() {
        return log10LikelihoodOfAFzero;
    }

    public double getLog10PosteriorOfAFzero() {
        return log10PosteriorOfAFzero;
    }

    public void reset() {
        log10MLE = log10MAP = log10LikelihoodOfAFzero = log10PosteriorOfAFzero = AlleleFrequencyCalculationModel.VALUE_NOT_CALCULATED;
        for ( int i = 0; i < alleleCountsOfMLE.length; i++ ) {
            alleleCountsOfMLE[i] = 0;
            alleleCountsOfMAP[i] = 0;
        }
        log10PosteriorMatrixValues.clear();
        log10PosteriorMatrixSum = null;
    }

    public void updateMLEifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        if ( log10LofK > log10MLE ) {
            log10MLE = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMLE[i] = alleleCountsForK[i];
        }
    }

    public void updateMAPifNeeded(final double log10LofK, final int[] alleleCountsForK) {
        log10PosteriorMatrixValues.add(log10LofK);
        if ( log10LofK > log10MAP ) {
            log10MAP = log10LofK;
            for ( int i = 0; i < alleleCountsForK.length; i++ )
                alleleCountsOfMAP[i] = alleleCountsForK[i];
        }
    }

    public void setLog10LikelihoodOfAFzero(final double log10LikelihoodOfAFzero) {
        this.log10LikelihoodOfAFzero = log10LikelihoodOfAFzero;
        if ( log10LikelihoodOfAFzero > log10MLE ) {
            log10MLE = log10LikelihoodOfAFzero;
            Arrays.fill(alleleCountsOfMLE, 0);
        }
    }

    public void setLog10PosteriorOfAFzero(final double log10PosteriorOfAFzero) {
        this.log10PosteriorOfAFzero = log10PosteriorOfAFzero;
        if ( log10PosteriorOfAFzero > log10MAP ) {
            log10MAP = log10PosteriorOfAFzero;
            Arrays.fill(alleleCountsOfMAP, 0);
        }
    }
}