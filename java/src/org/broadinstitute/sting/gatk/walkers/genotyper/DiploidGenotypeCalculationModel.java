/*
 * Copyright (c) 2010 The Broad Institute
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

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;

import java.util.*;

public class DiploidGenotypeCalculationModel extends JointEstimateGenotypeCalculationModel {
    /**
     * Should we enable the experimental genotype likelihood calculations?
     */
    protected boolean useExptGenotypeLikelihoods = false;

    protected DiploidGenotypeCalculationModel(boolean useExptGenotypeLikelihoods) {
        this.useExptGenotypeLikelihoods = useExptGenotypeLikelihoods;
    }

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();
    private HashMap<Character, AlleleFrequencyMatrix> AFMatrixMap = new HashMap<Character, AlleleFrequencyMatrix>();

    private enum GenotypeType { REF, HET, HOM }

    protected void initialize(char ref,
                              Map<String, StratifiedAlignmentContext> contexts,
                              StratifiedAlignmentContext.StratifiedContextType contextType) {
        // initialize the GenotypeLikelihoods
        GLs.clear();
        AFMatrixMap.clear();

        // for each alternate allele, create a new matrix
        for ( char alt : BaseUtils.BASES ) {
            if ( alt != ref )
                AFMatrixMap.put(alt, new AlleleFrequencyMatrix(contexts.size()));
        }

        // use flat priors for GLs
        DiploidGenotypePriors priors = new DiploidGenotypePriors();

        for ( String sample : contexts.keySet() ) {
            StratifiedAlignmentContext context = contexts.get(sample);
            ReadBackedPileup pileup = context.getContext(contextType).getBasePileup();

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = new GenotypeLikelihoods(UAC.baseModel, priors, UAC.defaultPlatform);

            GL.add(pileup, true, UAC.CAP_BASE_QUALITY);
            GLs.put(sample, GL);

            double[] posteriors = GL.getPosteriors();

            // for each alternate allele, fill the matrix
            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
            for ( char alt : BaseUtils.BASES ) {
                if ( alt != ref ) {
                    DiploidGenotype hetGenotype = DiploidGenotype.createDiploidGenotype(ref, alt);
                    DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(alt);
                    AFMatrixMap.get(alt).setLikelihoods(posteriors[refGenotype.ordinal()], posteriors[hetGenotype.ordinal()], posteriors[homGenotype.ordinal()], sample);
                }
            }
        }
    }

    protected void calculatelog10PofDgivenAFforAllF(char ref, char alt, int numFrequencies, Map<String, StratifiedAlignmentContext> contexts, StratifiedAlignmentContext.StratifiedContextType contextType) {

        AlleleFrequencyMatrix matrix = AFMatrixMap.get(alt);
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(alt);

        // first, calculate for AF=0 (no change to matrix)
        log10PofDgivenAFi[baseIndex][0] = matrix.getLikelihoodsOfFrequency();
        double maxLikelihoodSeen = log10PofDgivenAFi[baseIndex][0];
        int minAlleleFrequencyToTest = getMinAlleleFrequencyToTest();

        // for each minor allele frequency, calculate log10PofDgivenAFi
        for (int i = 1; i <= numFrequencies; i++) {
            // add one more alternatr allele
            matrix.incrementFrequency();

            // calculate new likelihoods
            log10PofDgivenAFi[baseIndex][i] = matrix.getLikelihoodsOfFrequency();

            // an optimization to speed up the calculation: if we are beyond the local maximum such
            //  that subsequent likelihoods won't factor into the confidence score, just quit
            if ( i >= minAlleleFrequencyToTest && maxLikelihoodSeen - log10PofDgivenAFi[baseIndex][i] > LOG10_OPTIMIZATION_EPSILON ) {
                ignoreAlleleFrequenciesAboveI(i, numFrequencies, baseIndex);
                return;
            }

            if ( log10PofDgivenAFi[baseIndex][i] > maxLikelihoodSeen )
                maxLikelihoodSeen = log10PofDgivenAFi[baseIndex][i];
        }
    }

    protected Map<String, Genotype> makeGenotypeCalls(char ref, char alt, int frequency, Map<String, StratifiedAlignmentContext> contexts, GenomeLoc loc) {
        HashMap<String, Genotype> calls = new HashMap<String, Genotype>();

        // set up some variables we'll need in the loop
        AlleleFrequencyMatrix matrix = AFMatrixMap.get(alt);
        Allele refAllele = new Allele(Character.toString(ref), true);
        Allele altAllele = new Allele(Character.toString(alt), false);
        DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
        DiploidGenotype hetGenotype = DiploidGenotype.createDiploidGenotype(ref, alt);
        DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(alt);

        for ( String sample : GLs.keySet() ) {

            // set the genotype and confidence
            Pair<Integer, Double> AFbasedGenotype = matrix.getGenotype(frequency, sample);

            ArrayList<Allele> myAlleles = new ArrayList<Allele>();
            if ( AFbasedGenotype.first == GenotypeType.REF.ordinal() ) {
                myAlleles.add(refAllele);
                myAlleles.add(refAllele);
            } else if ( AFbasedGenotype.first == GenotypeType.HET.ordinal() ) {
                myAlleles.add(refAllele);
                myAlleles.add(altAllele);
            } else { // ( AFbasedGenotype.first == GenotypeType.HOM.ordinal() )
                myAlleles.add(altAllele);
                myAlleles.add(altAllele);
            }

            CalledGenotype cg = new CalledGenotype(sample, myAlleles, AFbasedGenotype.second);
            cg.setLikelihoods(GLs.get(sample).getLikelihoods());
            cg.setReadBackedPileup(contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup());
            cg.putAttribute(VCFGenotypeRecord.DEPTH_KEY, contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size());

            double[] posteriors = GLs.get(sample).getPosteriors();
            cg.setPosteriors(posteriors);
            String GL = String.format("%.2f,%.2f,%.2f",
                    posteriors[refGenotype.ordinal()],
                    posteriors[hetGenotype.ordinal()],
                    posteriors[homGenotype.ordinal()]);
            cg.putAttribute(VCFGenotypeRecord.GENOTYPE_POSTERIORS_TRIPLET_KEY, GL);

            calls.put(sample, cg);
        }

        // output to beagle file if requested
        if ( beagleWriter != null ) {
            for ( String sample : samples ) {
                GenotypeLikelihoods gl = GLs.get(sample);
                if ( gl == null ) {
                    beagleWriter.print(" 0.0 0.0 0.0");
                    continue;
                }
                double[] likelihoods = gl.getLikelihoods();
                beagleWriter.print(' ');
                beagleWriter.print(String.format("%.6f", Math.pow(10, likelihoods[refGenotype.ordinal()])));
                beagleWriter.print(' ');
                beagleWriter.print(String.format("%.6f", Math.pow(10, likelihoods[hetGenotype.ordinal()])));
                beagleWriter.print(' ');
                beagleWriter.print(String.format("%.6f", Math.pow(10, likelihoods[homGenotype.ordinal()])));
            }
        }

        return calls;
    }


    protected class AlleleFrequencyMatrix {

        private double[][] matrix;    // allele frequency matrix
        private int[] indexes;        // matrix to maintain which genotype is active
        private int N;                // total frequencies
        private int frequency;        // current frequency

        // data structures necessary to maintain a list of the best genotypes and their scores
        private ArrayList<String> samples = new ArrayList<String>();
        private HashMap<Integer, HashMap<String, Pair<Integer, Double>>> samplesToGenotypesPerAF = new HashMap<Integer, HashMap<String, Pair<Integer, Double>>>();

        public AlleleFrequencyMatrix(int N) {
            this.N = N;
            frequency = 0;
            matrix = new double[N][3];
            indexes = new int[N];
            for (int i = 0; i < N; i++)
                indexes[i] = 0;
        }

        public void setLikelihoods(double AA, double AB, double BB, String sample) {
            int index = samples.size();
            samples.add(sample);
            matrix[index][GenotypeType.REF.ordinal()] = AA;
            matrix[index][GenotypeType.HET.ordinal()] = AB;
            matrix[index][GenotypeType.HOM.ordinal()] = BB;
        }

        public void incrementFrequency() {
            if ( frequency == 2 * N )
                throw new StingException("Frequency was incremented past N; how is this possible?");
            frequency++;

            double greedy = VALUE_NOT_CALCULATED;
            int greedyIndex = -1;
            for (int i = 0; i < N; i++) {

                if ( indexes[i] == GenotypeType.HET.ordinal() ) {
                    if ( matrix[i][GenotypeType.HOM.ordinal()] - matrix[i][GenotypeType.HET.ordinal()] > greedy ) {
                        greedy = matrix[i][GenotypeType.HOM.ordinal()] - matrix[i][GenotypeType.HET.ordinal()];
                        greedyIndex = i;
                    }
                }
                else if ( indexes[i] == GenotypeType.REF.ordinal() ) {
                    if ( matrix[i][GenotypeType.HET.ordinal()] - matrix[i][GenotypeType.REF.ordinal()] > greedy ) {
                        greedy = matrix[i][GenotypeType.HET.ordinal()] - matrix[i][GenotypeType.REF.ordinal()];
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
                throw new StingException("There is no best choice for a new alternate allele; how is this possible?");

            if ( indexes[greedyIndex] == GenotypeType.HET.ordinal() )
                indexes[greedyIndex] = GenotypeType.HOM.ordinal();
            else
                indexes[greedyIndex] = GenotypeType.HET.ordinal();
        }

        public double getLikelihoodsOfFrequency() {
            double likelihoods = 0.0;
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
                if ( genotype == GenotypeType.REF.ordinal() )
                    score = matrix[index][GenotypeType.REF.ordinal()] - Math.max(matrix[index][GenotypeType.HET.ordinal()], matrix[index][GenotypeType.HOM.ordinal()]);
                else if ( genotype == GenotypeType.HET.ordinal() )
                    score = matrix[index][GenotypeType.HET.ordinal()] - Math.max(matrix[index][GenotypeType.REF.ordinal()], matrix[index][GenotypeType.HOM.ordinal()]);
                else // ( genotype == GenotypeType.HOM.ordinal() )
                    score = matrix[index][GenotypeType.HOM.ordinal()] - Math.max(matrix[index][GenotypeType.REF.ordinal()], matrix[index][GenotypeType.HET.ordinal()]);

                samplesToGenotypes.put(sample, new Pair<Integer, Double>(genotype, Math.abs(score)));
                index++;
            }

            samplesToGenotypesPerAF.put(frequency, samplesToGenotypes);
        }
    }
}