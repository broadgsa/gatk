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
    private double[] log10Cache;

    protected ExactAFCalculationModel(int N, Logger logger, PrintStream verboseWriter) {
        super(N, logger, verboseWriter);

        log10Cache = new double[2*N];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k=1; k < 2*N; k++)
            log10Cache[k] = Math.log10(k);
    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, BiallelicGenotypeLikelihoods> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors,
                                int minFrequencyToCalculate) {


        int numSamples = GLs.size();
        int numChr = 2*numSamples;

        double[][] logYMatrix = new double[1+numSamples][1+numChr];

        for (int i=0; i <=numSamples; i++)
            for (int j=0; j <=numChr; j++)
                logYMatrix[i][j] = Double.NEGATIVE_INFINITY;

        //YMatrix[0][0] = 1.0;
        logYMatrix[0][0] = 0.0;
        int j=0;


        // call clear ref sites separately to speed computation
/*
        boolean isClearRefSite = true;

        int bestAFguess = 0;
        for ( String sample : GLs.keySet() ) {

            double[] genotypeLikelihoods = GLs.get(sample).getLikelihoods();
            if (!(genotypeLikelihoods[0] > genotypeLikelihoods[1] &&
                    genotypeLikelihoods[0] > genotypeLikelihoods[2]))
                isClearRefSite = false;
   
            bestAFguess += MathUtils.maxElementIndex(genotypeLikelihoods);

        }
     */
        for ( String sample : GLs.keySet() ) {
            j++;

            //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
            double[] genotypeLikelihoods = GLs.get(sample).getLikelihoods();
            //double logDenominator = Math.log10(2.0*j*(2.0*j-1));
            double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];

            // special treatment for k=0: iteration reduces to:
            //YMatrix[j][0] = YMatrix[j-1][0]*genotypeLikelihoods[GenotypeType.AA.ordinal()];
            logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[GenotypeType.AA.ordinal()];

            int k = 1;
            for (k=1; k <= 2*j; k++ ) {
   //             if (k > 3 && isClearRefSite)
   //                 break;


                //double num = (2.0*j-k)*(2.0*j-k-1)*YMatrix[j-1][k] * genotypeLikelihoods[GenotypeType.AA.ordinal()];
                double logNumerator[];
                logNumerator = new double[3];
                if (k < 2*j-1)
                    logNumerator[0] = log10Cache[2*j-k] + log10Cache[2*j-k-1] + logYMatrix[j-1][k] +
                        genotypeLikelihoods[GenotypeType.AA.ordinal()];
                else
                    logNumerator[0] = Double.NEGATIVE_INFINITY;


                if (k < 2*j)
                    logNumerator[1] = log10Cache[2*k] + log10Cache[2*j-k]+ logYMatrix[j-1][k-1] +
                        genotypeLikelihoods[GenotypeType.AB.ordinal()];
                else
                    logNumerator[1] = Double.NEGATIVE_INFINITY;

                if (k > 1)
                    logNumerator[2] = log10Cache[k] + log10Cache[k-1] + logYMatrix[j-1][k-2] +
                        genotypeLikelihoods[GenotypeType.BB.ordinal()];
                else
                    logNumerator[2] = Double.NEGATIVE_INFINITY;

                double logNum = softMax(logNumerator);

                //YMatrix[j][k] = num/den;
                logYMatrix[j][k] = logNum - logDenominator;
            }
      //      int kstart = k-1;
     //       for (k=kstart; k <= 2*j; k++)
       //         logYMatrix[j][k] = logYMatrix[j][kstart];


        }


        for (int k=0; k <= numChr; k++)
            log10AlleleFrequencyPosteriors[k] = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];


        // TODO: we really need to get rid of this and the minFrequencyToCalculate argument
        // it's possible that we need to calculate higher frequencies
        int maxAlleleFrequencyToTest = numChr;
        for (int i = maxAlleleFrequencyToTest; i <= minFrequencyToCalculate; i++)
            log10AlleleFrequencyPosteriors[i] = log10AlleleFrequencyPosteriors[maxAlleleFrequencyToTest];

    }


    double softMax(double[] vec) {
        // compute naively log10(10^x[0] + 10^x[1]+...)
        //return Math.log10(MathUtils.sumLog10(vec));

        //int maxInd = MathUtils.maxElementIndex(vec);
        //return vec[maxInd];

        // still more efficient: inline search for max. Hard-code fact that vec has length 3
        if (vec[0] > vec[1])
            if (vec[2] > vec[0])
                return vec[2];
            else
                return vec[0];
        else
            if (vec[2]>vec[1])
                return vec[2];
            else
                return vec[1];



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


            if (qual < 0) {
                // QUAL can be negative if the chosen genotype is not the most likely one individually.
                // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen on
                double[] normalized = MathUtils.normalizeFromLog10(posteriors);
                double chosenGenotype = normalized[bestGTguess];
                qual = -1.0 * Math.log10(1.0 - chosenGenotype);
            }

            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods(), UnifiedGenotyperV2.DEFAULT_GENOTYPE_LIKELIHOODS_KEY);
            attributes.put(likelihoods.getKey(), likelihoods.getAsString());

            calls.put(sample, new Genotype(sample, myAlleles, qual, null, attributes, false));

        }

        return calls;
    }

}
