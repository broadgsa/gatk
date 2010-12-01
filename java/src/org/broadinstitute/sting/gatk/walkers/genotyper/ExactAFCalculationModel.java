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

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.*;
import java.io.PrintStream;

public class ExactAFCalculationModel extends AlleleFrequencyCalculationModel {

    private boolean SIMPLE_GREEDY_GENOTYPER = false;
    private static final double[] log10Cache;
    private static final double[] jacobianLogTable;
    private static final int JACOBIAN_LOG_TABLE_SIZE = 101;
    private static final double JACOBIAN_LOG_TABLE_STEP = 0.1;
    private static final double MAX_JACOBIAN_TOLERANCE = 10.0;
    private static final int MAXN = 10000;

    static {
        log10Cache = new double[2*MAXN];
	    jacobianLogTable = new double[JACOBIAN_LOG_TABLE_SIZE];

        log10Cache[0] = Double.NEGATIVE_INFINITY;
        for (int k=1; k < 2*MAXN; k++)
            log10Cache[k] = Math.log10(k);

	for (int k=0; k < JACOBIAN_LOG_TABLE_SIZE; k++) {
	    jacobianLogTable[k] = Math.log10(1.0+Math.pow(10.0,-((double)k) 
						       * JACOBIAN_LOG_TABLE_STEP));

	}

    }

    protected ExactAFCalculationModel(int N, Logger logger, PrintStream verboseWriter) {
        super(N, logger, verboseWriter);

    }

    public void getLog10PNonRef(RefMetaDataTracker tracker,
                                ReferenceContext ref,
                                Map<String, Genotype> GLs,
                                double[] log10AlleleFrequencyPriors,
                                double[] log10AlleleFrequencyPosteriors) {


        int numSamples = GLs.size();
        int numChr = 2*numSamples;

        double[][] logYMatrix = new double[1+numSamples][1+numChr];

        for (int i=0; i <=numSamples; i++)
            for (int j=0; j <=numChr; j++)
                logYMatrix[i][j] = Double.NEGATIVE_INFINITY;

        //YMatrix[0][0] = 1.0;
        logYMatrix[0][0] = 0.0;
        int j=0;

        for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
            j++;

            if ( !sample.getValue().hasLikelihoods() )
                continue;

            //double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(GLs.get(sample).getLikelihoods());
            double[] genotypeLikelihoods = sample.getValue().getLikelihoods().getAsVector();
            //double logDenominator = Math.log10(2.0*j*(2.0*j-1));
            double logDenominator = log10Cache[2*j] + log10Cache[2*j-1];

            // special treatment for k=0: iteration reduces to:
            //YMatrix[j][0] = YMatrix[j-1][0]*genotypeLikelihoods[GenotypeType.AA.ordinal()];
            logYMatrix[j][0] = logYMatrix[j-1][0] + genotypeLikelihoods[GenotypeType.AA.ordinal()];

            for (int k=1; k <= 2*j; k++ ) {

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

        }


        for (int k=0; k <= numChr; k++)
            log10AlleleFrequencyPosteriors[k] = logYMatrix[j][k] + log10AlleleFrequencyPriors[k];

    }


    double softMax(double[] vec) {
        // compute naively log10(10^x[0] + 10^x[1]+...)
        //        return Math.log10(MathUtils.sumLog10(vec));

        // better approximation: do Jacobian logarithm function on data pairs
        double a = softMaxPair(vec[0],vec[1]);
        return softMaxPair(a,vec[2]);
    }

    double softMaxPair(double x, double y) {
        if (Double.isInfinite(x))
            return y;

        if (Double.isInfinite(y))
            return x;

        if (y >= x + MAX_JACOBIAN_TOLERANCE)
            return y;
        if (x >= y + MAX_JACOBIAN_TOLERANCE)
            return x;

        // OK, so |y-x| < tol: we use the following identity then:
        // we need to compute log10(10^x + 10^y)
        // By Jacobian logarithm identity, this is equal to
        // max(x,y) + log10(1+10^-abs(x-y))
        // we compute the second term as a table lookup
        // with integer quantization
        double diff = Math.abs(x-y);
        double t1 =x;
        if (y > x)
            t1 = y;
        // t has max(x,y)
        // we have pre-stored correction for 0,0.1,0.2,... 10.0
        int ind = (int)Math.round(diff/JACOBIAN_LOG_TABLE_STEP);
        double t2 = jacobianLogTable[ind];

        // gdebug+
        //double z =Math.log10(1+Math.pow(10.0,-diff));
        //System.out.format("x: %f, y:%f, app: %f, true: %f ind:%d\n",x,y,t2,z,ind);
        //gdebug-
        return t1+t2;
        // return Math.log10(Math.pow(10.0,x) + Math.pow(10.0,y));
    }



    /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    public Map<String, Genotype> assignGenotypes(VariantContext vc,
                                                 double[] log10AlleleFrequencyPosteriors,
                                                 int AFofMaxLikelihood) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        Allele refAllele = vc.getReference();
        Allele altAllele = vc.getAlternateAllele(0);

        Map<String, Genotype> GLs = vc.getGenotypes();
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

            for ( Map.Entry<String, Genotype> sample : GLs.entrySet() ) {
                if ( !sample.getValue().hasLikelihoods() )
                    continue;

                double[] likelihoods = sample.getValue().getLikelihoods().getAsVector();
                sampleIndices.add(sample.getKey());

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

        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();

        int startIdx = AFofMaxLikelihood;
        for (int k = sampleIdx; k > 0; k--) {
            int bestGTguess;
            String sample = sampleIndices.get(k-1);
            Genotype g = GLs.get(sample);
            if ( !g.hasLikelihoods() )
                continue;

            if (SIMPLE_GREEDY_GENOTYPER)
                bestGTguess = Utils.findIndexOfMaxEntry(g.getLikelihoods().getAsVector());
            else {
                int newIdx = tracebackArray[k][startIdx];
                bestGTguess = startIdx - newIdx;
                startIdx = newIdx;
            }

            ArrayList<Allele> myAlleles = new ArrayList<Allele>();

            double qual;
            double[] likelihoods = g.getLikelihoods().getAsVector();

            if (bestGTguess == 0) {
                myAlleles.add(refAllele);
                myAlleles.add(refAllele);
                qual = likelihoods[0] - Math.max(likelihoods[1], likelihoods[2]);
            } else if(bestGTguess == 1) {
                myAlleles.add(refAllele);
                myAlleles.add(altAllele);
                qual = likelihoods[1] - Math.max(likelihoods[0], likelihoods[2]);

            }  else {
                myAlleles.add(altAllele);
                myAlleles.add(altAllele);
                qual = likelihoods[2] - Math.max(likelihoods[1], likelihoods[0]);
            }


            if (qual < 0) {
                // QUAL can be negative if the chosen genotype is not the most likely one individually.
                // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen on
                double[] normalized = MathUtils.normalizeFromLog10(likelihoods);
                double chosenGenotype = normalized[bestGTguess];
                qual = -1.0 * Math.log10(1.0 - chosenGenotype);
            }

            calls.put(sample, new Genotype(sample, myAlleles, qual, null, g.getAttributes(), false));

        }

        return calls;
    }

}
