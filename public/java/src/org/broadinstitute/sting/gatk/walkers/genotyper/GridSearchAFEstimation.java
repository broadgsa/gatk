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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
import java.util.*;

public class GridSearchAFEstimation extends AlleleFrequencyCalculationModel {

    // for use in optimizing the P(D|AF) calculations:
    // how much off from the max likelihoods do we need to be before we can quit calculating?
    protected static final double LOG10_OPTIMIZATION_EPSILON = 8.0;

    private AlleleFrequencyMatrix AFMatrix;

    protected GridSearchAFEstimation(UnifiedArgumentCollection UAC, int N, Logger logger, PrintStream verboseWriter) {
        super(UAC, N, logger, verboseWriter);
        AFMatrix = new AlleleFrequencyMatrix(N);
    }

    protected void getLog10PNonRef(Map<String, Genotype> GLs, List<Allele> alleles,
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
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
    protected Map<String, Genotype> assignGenotypes(VariantContext vc,
                                                    double[] log10AlleleFrequencyPosteriors,
                                                    int AFofMaxLikelihood) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        Allele refAllele = vc.getReference();
        Allele altAllele = vc.getAlternateAllele(0);
        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();

        // first, the potential alt calls
        for ( String sample : AFMatrix.getSamples() ) {
            Genotype g = vc.getGenotype(sample);

            // set the genotype and confidence
            Pair<Integer, Double> AFbasedGenotype = AFMatrix.getGenotype(AFofMaxLikelihood, sample);
            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            if ( AFbasedGenotype.first == GenotypeType.AA.ordinal() ) {
                myAlleles.add(refAllele);
                myAlleles.add(refAllele);
            } else if ( AFbasedGenotype.first == GenotypeType.AB.ordinal() ) {
                myAlleles.add(refAllele);
                myAlleles.add(altAllele);
            } else { // ( AFbasedGenotype.first == GenotypeType.BB.ordinal() )
                myAlleles.add(altAllele);
                myAlleles.add(altAllele);
            }

            calls.put(sample, new Genotype(sample, myAlleles, AFbasedGenotype.second, null, g.getAttributes(), false));
        }

        return calls;
    }

    private void initializeAFMatrix(Map<String, Genotype> GLs) {
        AFMatrix.clear();

        for ( Genotype g : GLs.values() ) {
            if ( g.hasLikelihoods() )
                AFMatrix.setLikelihoods(g.getLikelihoods().getAsVector(), g.getSampleName());
        }
    }

    protected static class AlleleFrequencyMatrix {

        private double[][] matrix;    // allele frequency matrix
        private int[] indexes;        // matrix to maintain which genotype is active
        private int maxN;             // total possible frequencies in data
        private int frequency;        // current frequency

        // data structures necessary to maintain a list of the best genotypes and their scores
        private ArrayList<String> samples = new ArrayList<String>();
        private HashMap<Integer, HashMap<String, Pair<Integer, Double>>> samplesToGenotypesPerAF = new HashMap<Integer, HashMap<String, Pair<Integer, Double>>>();

        public AlleleFrequencyMatrix(int N) {
            maxN = N;
            matrix = new double[N][3];
            indexes = new int[N];
            clear();
        }

        public List<String> getSamples() { return samples; }

        public void clear() {
            frequency = 0;
            for (int i = 0; i < maxN; i++)
                indexes[i] = 0;
            samples.clear();
            samplesToGenotypesPerAF.clear();
        }

        public void setLikelihoods(double[] GLs, String sample) {
            int index = samples.size();
            samples.add(sample);
            matrix[index][GenotypeType.AA.ordinal()] = GLs[0];
            matrix[index][GenotypeType.AB.ordinal()] = GLs[1];
            matrix[index][GenotypeType.BB.ordinal()] = GLs[2];
        }

        public void incrementFrequency() {
            int N = samples.size();
            if ( frequency == 2 * N )
                throw new ReviewedStingException("Frequency was incremented past N; how is this possible?");
            frequency++;

            double greedy = VALUE_NOT_CALCULATED;
            int greedyIndex = -1;
            for (int i = 0; i < N; i++) {

                if ( indexes[i] == GenotypeType.AB.ordinal() ) {
                    if ( matrix[i][GenotypeType.BB.ordinal()] - matrix[i][GenotypeType.AB.ordinal()] > greedy ) {
                        greedy = matrix[i][GenotypeType.BB.ordinal()] - matrix[i][GenotypeType.AB.ordinal()];
                        greedyIndex = i;
                    }
                }
                else if ( indexes[i] == GenotypeType.AA.ordinal() ) {
                    if ( matrix[i][GenotypeType.AB.ordinal()] - matrix[i][GenotypeType.AA.ordinal()] > greedy ) {
                        greedy = matrix[i][GenotypeType.AB.ordinal()] - matrix[i][GenotypeType.AA.ordinal()];
                        greedyIndex = i;
                    }
                    // note that we currently don't bother with breaking ties between samples
                    // (which would be done by looking at the HOM_VAR value) because it's highly
                    // unlikely that a collision will both occur and that the difference will
                    // be significant at HOM_VAR...
                }
                // if this person is already hom var, he can't add another alternate allele
                // so we can ignore that case
            }
            if ( greedyIndex == -1 )
                throw new ReviewedStingException("There is no best choice for a new alternate allele; how is this possible?");

            if ( indexes[greedyIndex] == GenotypeType.AB.ordinal() )
                indexes[greedyIndex] = GenotypeType.BB.ordinal();
            else
                indexes[greedyIndex] = GenotypeType.AB.ordinal();
        }

        public double getLikelihoodsOfFrequency() {
            double likelihoods = 0.0;
            int N = samples.size();
            for (int i = 0; i < N; i++)
                likelihoods += matrix[i][indexes[i]];

            /*
            System.out.println(frequency);
            for (int i = 0; i < N; i++) {
                System.out.print(samples.get(i));
                for (int j=0; j < 3; j++) {
                    System.out.print(String.valueOf(matrix[i][j]));
                    System.out.print(indexes[i] == j ? "* " : " ");
                }
                System.out.println();
            }
            System.out.println(likelihoods);
            System.out.println();
            */

            recordGenotypes();

            return likelihoods;
        }

        public Pair<Integer, Double> getGenotype(int frequency, String sample) {
            return samplesToGenotypesPerAF.get(frequency).get(sample);
        }

        private void recordGenotypes() {
            HashMap<String, Pair<Integer, Double>> samplesToGenotypes = new HashMap<String, Pair<Integer, Double>>();

            int index = 0;
            for ( String sample : samples ) {
                int genotype = indexes[index];

                double score;

                int maxEntry = MathUtils.maxElementIndex(matrix[index]);
                // if the max value is for the most likely genotype, we can compute next vs. next best
                if ( genotype == maxEntry ) {
                    if ( genotype == GenotypeType.AA.ordinal() )
                        score = matrix[index][genotype] - Math.max(matrix[index][GenotypeType.AB.ordinal()], matrix[index][GenotypeType.BB.ordinal()]);
                    else if ( genotype == GenotypeType.AB.ordinal() )
                        score = matrix[index][genotype] - Math.max(matrix[index][GenotypeType.AA.ordinal()], matrix[index][GenotypeType.BB.ordinal()]);
                    else // ( genotype == GenotypeType.HOM.ordinal() )
                        score = matrix[index][genotype] - Math.max(matrix[index][GenotypeType.AA.ordinal()], matrix[index][GenotypeType.AB.ordinal()]);
                }
                // otherwise, we need to calculate the probability of the genotype
                else {
                    double[] normalized = MathUtils.normalizeFromLog10(matrix[index]);
                    double chosenGenotype = normalized[genotype];
                    score = -1.0 * Math.log10(1.0 - chosenGenotype);
                }

                samplesToGenotypes.put(sample, new Pair<Integer, Double>(genotype, Math.abs(score)));
                index++;
            }

            samplesToGenotypesPerAF.put(frequency, samplesToGenotypes);
        }
    }
}
