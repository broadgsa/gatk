package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

public class DiploidGenotypeCalculationModel extends JointEstimateGenotypeCalculationModel {

    protected DiploidGenotypeCalculationModel() {}

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();
    private HashMap<Character, AlleleFrequencyMatrix> AFMatrixMap = new HashMap<Character, AlleleFrequencyMatrix>();

    private enum GenotypeType { REF, HET, HOM }

    protected void initialize(char ref, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
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

        int index = 0;
        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = context.getContext(contextType).getPileup();

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = new GenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.add(pileup, true);
            GLs.put(sample, GL);

            double[] posteriors = GL.getPosteriors();

            // for each alternate allele, fill the matrix
            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
            for ( char alt : BaseUtils.BASES ) {
                if ( alt != ref ) {
                    DiploidGenotype hetGenotype = ref < alt ? DiploidGenotype.valueOf(String.valueOf(ref) + String.valueOf(alt)) : DiploidGenotype.valueOf(String.valueOf(alt) + String.valueOf(ref));
                    DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(alt);
                    AFMatrixMap.get(alt).setLikelihoods(posteriors[refGenotype.ordinal()], posteriors[hetGenotype.ordinal()], posteriors[homGenotype.ordinal()], index);
                }
            }
            index++;
        }
    }

    protected void calculatelog10PofDgivenAFforAllF(char ref, char alt, int numFrequencies, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {

        AlleleFrequencyMatrix matrix = AFMatrixMap.get(alt);
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(alt);

        // first, calculate for AF=0 (no change to matrix)
        log10PofDgivenAFi[baseIndex][0] = matrix.getLikelihoodsOfFrequency();

        // for each minor allele frequency, calculate log10PofDgivenAFi
        for (int i = 1; i <= numFrequencies; i++) {
            // add one more alternatr allele
            matrix.incrementFrequency();

            // calculate new likelihoods
            log10PofDgivenAFi[baseIndex][i] = matrix.getLikelihoodsOfFrequency();
        }
    }

    protected List<Genotype> makeGenotypeCalls(char ref, char alt, HashMap<String, AlignmentContextBySample> contexts, GenomeLoc loc) {
        ArrayList<Genotype> calls = new ArrayList<Genotype>();

        for ( String sample : GLs.keySet() ) {

            // create the call
            GenotypeCall call = GenotypeWriterFactory.createSupportedGenotypeCall(OUTPUT_FORMAT, ref, loc);

            if ( call instanceof ReadBacked ) {
                ReadBackedPileup pileup = contexts.get(sample).getContext(StratifiedContext.OVERALL).getPileup();
                ((ReadBacked)call).setPileup(pileup);
            }
            if ( call instanceof SampleBacked ) {
                ((SampleBacked)call).setSampleName(sample);
            }
            if ( call instanceof LikelihoodsBacked ) {
                ((LikelihoodsBacked)call).setLikelihoods(GLs.get(sample).getLikelihoods());
            }
            if ( call instanceof PosteriorsBacked ) {
                ((PosteriorsBacked)call).setPosteriors(GLs.get(sample).getPosteriors());
            }
            if ( call instanceof AlleleConstrainedGenotype ) {
                ((AlleleConstrainedGenotype)call).setAlternateAllele(alt);
            }

            calls.add(call);
        }

        return calls;
    }

    protected class AlleleFrequencyMatrix {

        private double[][] matrix;
        private int[] indexes;
        private int N;
        private int frequency;

        public AlleleFrequencyMatrix(int N) {
            this.N = N;
            frequency = 0;
            matrix = new double[N][3];
            indexes = new int[N];
            for (int i = 0; i < N; i++)
                indexes[i] = 0;
        }

        public void setLikelihoods(double AA, double AB, double BB, int index) {
            matrix[index][GenotypeType.REF.ordinal()] = AA;
            matrix[index][GenotypeType.HET.ordinal()] = AB;
            matrix[index][GenotypeType.HOM.ordinal()] = BB;
        }

        public void incrementFrequency() {
            if ( frequency == 2 * N )
                throw new StingException("Frequency was incremented past N; how is this possible?");
            frequency++;

            double greedy = -1.0 * Double.MAX_VALUE;
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

            //verboseWriter.write(frequency + "\n");
            //for (int i = 0; i < N; i++) {
            //    for (int j=0; j < 3; j++) {
            //        verboseWriter.write(String.valueOf(matrix[i][j]));
            //        verboseWriter.write(indexes[i] == j ? "* " : " ");
            //    }
            //    verboseWriter.write("\n");
            //}
            //verboseWriter.write(likelihoods + "\n\n");

            return likelihoods;
        }
    }
}