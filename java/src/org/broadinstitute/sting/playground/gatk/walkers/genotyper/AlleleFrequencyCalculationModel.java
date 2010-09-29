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
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.vcf.VCFConstants;

import java.io.PrintStream;
import java.util.*;


/**
 * The model representing how we calculate a genotype given the priors and a pile
 * of bases and quality scores
 */
public abstract class AlleleFrequencyCalculationModel implements Cloneable {

    public enum Model {
        EXACT,
        GRID_SEARCH
    }

    protected int N;
    protected AlleleFrequencyMatrix AFMatrix;
    protected Set<BiallelicGenotypeLikelihoods> refCalls;
    
    protected Logger logger;
    protected PrintStream verboseWriter;

    private int minAlleleFrequencyToTest;

    protected AlleleFrequencyCalculationModel(int N, Logger logger, PrintStream verboseWriter) {
        this.N = N;
        this.logger = logger;
        this.verboseWriter = verboseWriter;
        AFMatrix = new AlleleFrequencyMatrix(N);
        refCalls = new HashSet<BiallelicGenotypeLikelihoods>();
    }

    /**
     * Must be overridden by concrete subclasses
     * @param tracker                         rod data
     * @param ref                             reference context
     * @param GLs                             genotype likelihoods
     * @param log10AlleleFrequencyPriors      priors
     * @param log10AlleleFrequencyPosteriors  array (pre-allocated) to store results
     */
    public abstract void getLog10PNonRef(RefMetaDataTracker tracker,
                                         ReferenceContext ref,
                                         Map<String, BiallelicGenotypeLikelihoods> GLs,
                                         double[] log10AlleleFrequencyPriors,
                                         double[] log10AlleleFrequencyPosteriors);

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
        initializeAFMatrix(GLs);

        // increment the grid
        for (int i = 1; i <= AFofMaxLikelihood; i++) {
            // add one more alternate allele
            AFMatrix.incrementFrequency();
        }

        return generateCalls(contexts, GLs, AFofMaxLikelihood);
    }

    // TODO: get rid of this optimization, it is wrong!
    protected int getMinAlleleFrequencyToTest() {
        return minAlleleFrequencyToTest;
    }

    protected void setMinAlleleFrequencyToTest(int minAF) {
        minAlleleFrequencyToTest = minAF;
    }

    protected Map<String, Genotype> generateCalls(Map<String, StratifiedAlignmentContext> contexts,
                                                  Map<String, BiallelicGenotypeLikelihoods> GLs,
                                                  int frequency) {
        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();

        // first, the potential alt calls
        for ( String sample : AFMatrix.getSamples() ) {
            BiallelicGenotypeLikelihoods GL = GLs.get(sample);
            Allele alleleA = GL.getAlleleA();
            Allele alleleB = GL.getAlleleB();

            // set the genotype and confidence
            Pair<Integer, Double> AFbasedGenotype = AFMatrix.getGenotype(frequency, sample);
            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            if ( AFbasedGenotype.first == GenotypeType.AA.ordinal() ) {
                myAlleles.add(alleleA);
                myAlleles.add(alleleA);
            } else if ( AFbasedGenotype.first == GenotypeType.AB.ordinal() ) {
                myAlleles.add(alleleA);
                myAlleles.add(alleleB);
            } else { // ( AFbasedGenotype.first == GenotypeType.BB.ordinal() )
                myAlleles.add(alleleB);
                myAlleles.add(alleleB);
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.DEPTH_KEY, contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size());

            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods.getAsString());

            calls.put(sample, new Genotype(sample, myAlleles, AFbasedGenotype.second, null, attributes, false));
        }

        // now, the clearly ref calls
        for ( BiallelicGenotypeLikelihoods GL : refCalls ) {
            String sample = GL.getSample();

            Allele ref = GL.getAlleleA();
            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            myAlleles.add(ref);
            myAlleles.add(ref);

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            attributes.put(VCFConstants.DEPTH_KEY, contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size());

            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods.getAsString());

            double GQ = GL.getAALikelihoods() - Math.max(GL.getABLikelihoods(), GL.getBBLikelihoods());

            calls.put(sample, new Genotype(sample, myAlleles, GQ, null, attributes, false));

        }

        return calls;
    }

    protected void initializeAFMatrix(Map<String, BiallelicGenotypeLikelihoods> GLs) {
        refCalls.clear();
        AFMatrix.clear();

        for ( BiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( isClearRefCall(GL) ) {
                refCalls.add(GL);
            } else {
                AFMatrix.setLikelihoods(GL.getPosteriors(), GL.getSample());
            }
        }
    }

    private boolean isClearRefCall(BiallelicGenotypeLikelihoods GL) {
        if ( GL.getAlleleA().isNonReference() )
            return false;

        double[] likelihoods = GL.getLikelihoods();
        return ( likelihoods[0] > likelihoods[1] && likelihoods[0] > likelihoods[2]);
    }

    protected class CalculatedAlleleFrequency {

        public double log10PNonRef;
        public int altAlleles;

        public CalculatedAlleleFrequency(double log10PNonRef, int altAlleles) {
            this.log10PNonRef = log10PNonRef;
            this.altAlleles = altAlleles;
        }
    }

    private enum GenotypeType { AA, AB, BB }
    protected static final double VALUE_NOT_CALCULATED = -1.0 * Double.MAX_VALUE;

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

        public void setLikelihoods(double AA, double AB, double BB, String sample) {
            int index = samples.size();
            samples.add(sample);
            matrix[index][GenotypeType.AA.ordinal()] = AA;
            matrix[index][GenotypeType.AB.ordinal()] = AB;
            matrix[index][GenotypeType.BB.ordinal()] = BB;
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