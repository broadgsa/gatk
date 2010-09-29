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
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.*;
import java.io.PrintStream;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {

    private int maxNumIterations = 50;
    private double tol = 1e-4;
    private boolean DEBUGOUT = true;

    protected ExactAFCalculationModel(int N, Logger logger, PrintStream verboseWriter) {
        super(N, logger, verboseWriter);
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, BiallelicGenotypeLikelihoods> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors) {

        // Math requires linear math to make efficient updates.
        double[] alleleFrequencyPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPriors);

        // now that we have genotype likelihoods for each sample, we can start refining allele frequency estimate
        double meanAF = computeMeanAF(alleleFrequencyPosteriors);
        if (DEBUGOUT)
            System.out.format("Initial Mean AF: %5.6f\n", meanAF);


        for (int numIterations=0; numIterations < maxNumIterations; numIterations++) {
            double oldMeanAF = meanAF;
            alleleFrequencyPosteriors = updateAFEstimate(GLs, alleleFrequencyPosteriors);
            meanAF = computeMeanAF(alleleFrequencyPosteriors);

            if (DEBUGOUT)
                System.out.format("Mean AF: %5.4f. PVariant: %5.5f\n", meanAF,1.0-alleleFrequencyPosteriors[0]);

            if (Math.abs(meanAF-oldMeanAF) < tol)
                break;

        }

        for (int k=0; k < alleleFrequencyPosteriors.length; k++)
            log10AlleleFrequencyPosteriors[k] = Math.log10(alleleFrequencyPosteriors[k]);

    }

    private double[] updateAFEstimate(Map<String, BiallelicGenotypeLikelihoods> GLs, double[] alleleFrequencyPriors) {

        int numSamples = GLs.size();
        int numChr = 2*numSamples;

        double[][] YMatrix = new double[1+numSamples][1+numChr];
        YMatrix[0][0] = 1.0;
        int j=0;

        for ( String sample : GLs.keySet() ) {
            j++;

            double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
            double sum = 0.0;
            double den = 2.0*j*(2.0*j-1);

            for (int k=0; k <= 2*j; k++ ) {
                double tmp = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[2];
                if (k > 0)
                    tmp += (2.0*k)*(2.0*j-k)*YMatrix[j-1][k-1]*genotypeLikelihoods[1];
                if (k > 1)
                    tmp += k*(k-1)*YMatrix[j-1][k-2]*genotypeLikelihoods[0];

                YMatrix[j][k] = tmp/den;

                sum += YMatrix[j][k];
            }
            // renormalize row
            for (int k=0; k <= 2*j; k++ )
                YMatrix[j][k] /= sum;

        }
        double sum = 0.0;
        double[] newAF = new double[alleleFrequencyPriors.length];

        for (int k=0; k <= numChr; k++) {
            double prod = YMatrix[j][numChr-k] * alleleFrequencyPriors[k];
            newAF[k] = prod;
            sum += prod;
        }
        //renormalize now
        for (int k=0; k < newAF.length; k++)
            newAF[k] /= sum;

        return newAF;


    }

    private double computeMeanAF(double[] afVector) {
        // get now new site AF estimate
        double sum = 0.0;
        for (int k=0; k < afVector.length; k++)
            sum += (double)k * afVector[k];

        return  sum/afVector.length;

    }

}