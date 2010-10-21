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
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.*;
import java.io.PrintStream;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {

    private boolean DEBUGOUT = false;
    private boolean SIMPLE_GREEDY_GENOTYPER = false;
    private static final double EPS = 1e-300;
    private static final double LOGEPS = -300;

    protected ExactAFCalculationModel(int N, Logger logger, PrintStream verboseWriter) {
        super(N, logger, verboseWriter);
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, BiallelicGenotypeLikelihoods> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors,
                                int minFrequencyToCalculate) {

        // Math requires linear math to make efficient updates.
        double[] alleleFrequencyPriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPriors);
        double[] alleleFrequencyPosteriors;
 

        alleleFrequencyPosteriors = updateAFEstimate(GLs, alleleFrequencyPriors);

        if (DEBUGOUT) {
            double meanAF = computeMeanAF(alleleFrequencyPosteriors);
            System.out.format("Mean AF: %5.4f. PVariant: %5.5f\n", meanAF,1.0-alleleFrequencyPosteriors[0]);
        }


        for (int k=0; k < alleleFrequencyPosteriors.length; k++) {
            if (alleleFrequencyPosteriors[k] > 1-EPS)
                log10AlleleFrequencyPosteriors[k] = -EPS;
            else if (alleleFrequencyPosteriors[k] < EPS)
                log10AlleleFrequencyPosteriors[k] = LOGEPS;
            else
                log10AlleleFrequencyPosteriors[k] = Math.log10(alleleFrequencyPosteriors[k]);
        }
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
                double tmp = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[GenotypeType.AA.ordinal()];
                if (k > 0)
                    tmp += (2.0*k)*(2.0*j-k)*YMatrix[j-1][k-1]*genotypeLikelihoods[GenotypeType.AB.ordinal()];
                if (k > 1)
                    tmp += k*(k-1)*YMatrix[j-1][k-2]*genotypeLikelihoods[GenotypeType.BB.ordinal()];

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
            double prod = YMatrix[j][k] * alleleFrequencyPriors[k];
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

    /**
     * Can be overridden by concrete subclasses
     * @param contexts             alignment contexts
     * @param GLs                  genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    public Map<String, Genotype> assignGenotypes(Map<String, StratifiedAlignmentContext> contexts,
                                                 Map<String, BiallelicGenotypeLikelihoods> GLs,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int AFofMaxLikelihood) {
        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();


        double[][] pathMetricArray = new double[GLs.size()+1][AFofMaxLikelihood+1];
        int[][] tracebackArray = new int[GLs.size()+1][AFofMaxLikelihood+1];

        ArrayList<String> sampleIndices = new ArrayList<String>();
        int sampleIdx = 0;

        // todo - optimize initialization
        for (int k=0; k <= AFofMaxLikelihood; k++)
            for (int j=0; j <= GLs.size(); j++)
                pathMetricArray[j][k] = -1e30;

        pathMetricArray[0][0] = 0.0;

        if (SIMPLE_GREEDY_GENOTYPER) {
            sampleIndices.addAll(GLs.keySet());
            sampleIdx = GLs.size();
        }
        else {

            for (String sample: GLs.keySet()) {
                sampleIndices.add(sample);

                double[] likelihoods = GLs.get(sample).getLikelihoods();

                for (int k=0; k <= AFofMaxLikelihood; k++) {

                    double bestMetric = pathMetricArray[sampleIdx][k] + likelihoods[0];
                    int bestIndex = k;

                    if (k>0) {
                        double m2 =  pathMetricArray[sampleIdx][k-1] + likelihoods[1];
                        if (m2 > bestMetric) {
                            bestMetric = m2;
                            bestIndex  = k-1;
                        }
                    }

                    if (k>1) {
                        double m2 =  pathMetricArray[sampleIdx][k-2] + likelihoods[2];
                        if (m2 > bestMetric) {
                            bestMetric = m2;
                            bestIndex  = k-2;
                        }
                    }

                    pathMetricArray[sampleIdx+1][k] = bestMetric;
                    tracebackArray[sampleIdx+1][k] = bestIndex;
                }
                sampleIdx++;
            }
        }

        int startIdx = AFofMaxLikelihood;
        for (int k = sampleIdx; k > 0; k--) {
            int bestGTguess;
            String sample = sampleIndices.get(k-1);
            BiallelicGenotypeLikelihoods GL = GLs.get(sample);
            Allele alleleA = GL.getAlleleA();
            Allele alleleB = GL.getAlleleB();

            if (SIMPLE_GREEDY_GENOTYPER)
                bestGTguess = Utils.findIndexOfMaxEntry(GLs.get(sample).getLikelihoods());
            else {
                int newIdx = tracebackArray[k][startIdx];
                bestGTguess = startIdx - newIdx;
                startIdx = newIdx;
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            AlignmentContext context = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

            if (context.hasBasePileup())
                attributes.put(VCFConstants.DEPTH_KEY, getFilteredDepth(context.getBasePileup()));

            else if (context.hasExtendedEventPileup())
                attributes.put(VCFConstants.DEPTH_KEY, getFilteredDepth(context.getExtendedEventPileup()));

            double qual;
            double[] posteriors = GLs.get(sample).getPosteriors();

            if (bestGTguess == 0) {
                myAlleles.add(alleleA);
                myAlleles.add(alleleA);
                qual = posteriors[0] - Math.max(posteriors[1],posteriors[2]);
            } else if(bestGTguess == 1) {
                myAlleles.add(alleleA);
                myAlleles.add(alleleB);
                qual = posteriors[1] - Math.max(posteriors[0],posteriors[2]);

            }  else {
                myAlleles.add(alleleB);
                myAlleles.add(alleleB);
                qual = posteriors[2] - Math.max(posteriors[1],posteriors[0]);
            }

            if (qual <= 0.0) {
                qual = 0.0;
             }

            attributes.put(VCFConstants.GENOTYPE_QUALITY_KEY,String.format("%4.2f", 10*qual));

            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods.getAsString());
            calls.put(sample, new Genotype(sample, myAlleles, qual, null, attributes, false));

        }

        return calls;
    }

}
