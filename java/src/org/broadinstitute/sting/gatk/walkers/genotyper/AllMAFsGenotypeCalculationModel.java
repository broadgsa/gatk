package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

public class AllMAFsGenotypeCalculationModel extends GenotypeCalculationModel {

    protected AllMAFsGenotypeCalculationModel() {}

    // because the null allele frequencies are constant for a given N,
    // we cache the results to avoid having to recompute everything
    private HashMap<Integer, double[]> nullAlleleFrequencyCache = new HashMap<Integer, double[]>();

    // because the Hardy-Weinberg values for a given frequency are constant,
    // we cache the results to avoid having to recompute everything
    private HashMap<Double, double[]> hardyWeinbergValueCache = new HashMap<Double, double[]>();

    // the allele frequency priors
    private double[] alleleFrequencyPriors;

    // the allele frequency posteriors and P(f>0) mapped by alternate allele
    private HashMap<Character, double[]> alleleFrequencyPosteriors = new HashMap<Character, double[]>();
    private HashMap<Character, Double> PofFs = new HashMap<Character, Double>();

    // the minimum and actual number of points in our allele frequency estimation
    private static final int MIN_ESTIMATION_POINTS = 20; //1000;
    private int frequencyEstimationPoints;

    // the GenotypeLikelihoods map
    private HashMap<String, AlleleSpecificGenotypeLikelihoods> GLs = new HashMap<String, AlleleSpecificGenotypeLikelihoods>();


    public Pair<List<GenotypeCall>, GenotypeMetaData> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the context for each sample, overall and separated by strand
        HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context);
        if ( contexts == null )
            return null;

        calculate(ref, contexts, priors, StratifiedContext.OVERALL);


        //double lod = overall.getPofD() - overall.getPofNull();
        //logger.debug("lod=" + lod);

        // calculate strand score
        //EMOutput forward = calculate(ref, contexts, priors, StratifiedContext.FORWARD);
        //EMOutput reverse = calculate(ref, contexts, priors, StratifiedContext.REVERSE);
        //double forwardLod = (forward.getPofD() + reverse.getPofNull()) - overall.getPofNull();
        //double reverseLod = (reverse.getPofD() + forward.getPofNull()) - overall.getPofNull();
        //logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);
        //double strandScore = Math.max(forwardLod - lod, reverseLod - lod);

        //logger.debug(String.format("LOD=%f, SLOD=%f", lod, strandScore));

        // generate the calls
        //GenotypeMetaData metadata = new GenotypeMetaData(lod, strandScore, overall.getMAF());
        //return new Pair<List<GenotypeCall>, GenotypeMetaData>(genotypeCallsFromGenotypeLikelihoods(overall, ref, contexts), metadata);
        return null;
    }

    private void calculate(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {

        initializeAlleleFrequencies(contexts.size());

        initializeGenotypeLikelihoods(ref, contexts, priors, contextType);

        calculateAlleleFrequencyPosteriors(ref);
    }

    private void initializeAlleleFrequencies(int numSamples) {

        // calculate the number of estimation points to use:
        // it's either MIN_ESTIMATION_POINTS or 2N if that's larger
        // (add 1 for allele frequency of zero)
        frequencyEstimationPoints = Math.max(MIN_ESTIMATION_POINTS, 2 * numSamples) + 1;

        // set up the allele frequency priors
        alleleFrequencyPriors = getNullalleleFrequencyPriors(frequencyEstimationPoints);

        for (int i = 1; i < frequencyEstimationPoints; i++)
            logger.debug("Null allele frequency for MAF=" + i + ": " + alleleFrequencyPriors[i]);
    }

    private double[] getNullalleleFrequencyPriors(int N) {
        double[] AFs = nullAlleleFrequencyCache.get(N);

        // if it hasn't been calculated yet, do so now
        if ( AFs == null ) {

            // calculate sum(1/i)
            double denominator = 0.0;
            for (int i = 1; i < N; i++)
                denominator += 1.0 / (double)i;

            // set up delta
            double delta = 1.0 / denominator;

            // calculate the null allele frequencies
            AFs = new double[N];
            for (int i = 1; i < N; i++)
                AFs[i] = Math.log10(delta / (double)i);

            nullAlleleFrequencyCache.put(N, AFs);
        }

        return AFs.clone();
    }

    private void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, StratifiedContext contextType) {
        GLs.clear();

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = GenotypeLikelihoodsFactory.makeGenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.setVerbose(VERBOSE);
            GL.add(pileup, true);

            GLs.put(sample, new AlleleSpecificGenotypeLikelihoods(ref, GL));
        }
    }

    private void calculateAlleleFrequencyPosteriors(char ref) {

        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;

            double[] AFPosteriors = new double[frequencyEstimationPoints];

            for ( AlleleSpecificGenotypeLikelihoods GL : GLs.values() ) {
                double[] PofDgivenG = GL.getNormalizedPosteriors(altAllele);

                // calculate the posterior weighted frequencies for this sample
                for (int i = 1; i < frequencyEstimationPoints; i++) {
                    double f = (double)i / (double)(frequencyEstimationPoints-1);
                    double[] hardyWeinberg = getHardyWeinbergValues(f);

                    double PofDgivenAFi = 0.0;
                    PofDgivenAFi += hardyWeinberg[GenotypeType.REF.ordinal()] * PofDgivenG[GenotypeType.REF.ordinal()];
                    PofDgivenAFi += hardyWeinberg[GenotypeType.HET.ordinal()] * PofDgivenG[GenotypeType.HET.ordinal()];
                    PofDgivenAFi += hardyWeinberg[GenotypeType.HOM.ordinal()] * PofDgivenG[GenotypeType.HOM.ordinal()];

                    AFPosteriors[i] += Math.log10(PofDgivenAFi);
                }
            }

            for (int i = 1; i < frequencyEstimationPoints; i++)
                logger.debug("P(D|AF=" + i + ") for alt allele " + altAllele + ": " + AFPosteriors[i]);

            // 1) multiply by null allele frequency priors to get AF posteriors
            // 2) for precision purposes, we need to subtract the largest posterior value
            // 3) convert back to normal space
            double maxValue = findMaxEntry(AFPosteriors, true);
            for (int i = 1; i < frequencyEstimationPoints; i++)
                AFPosteriors[i] = Math.pow(10, alleleFrequencyPriors[i] + AFPosteriors[i] - maxValue);

            // normalize
            double sum = 0.0;
            for (int i = 1; i < frequencyEstimationPoints; i++)
                sum += AFPosteriors[i];
            for (int i = 1; i < frequencyEstimationPoints; i++)
                AFPosteriors[i] /= sum;

            for (int i = 1; i < frequencyEstimationPoints; i++)
                logger.debug("Allele frequency posterior estimate for alt allele " + altAllele + " MAF=" + i + ": " + AFPosteriors[i]);
            logger.debug("P(f>0) for alt allele " + altAllele + ": " + sum);

            alleleFrequencyPosteriors.put(altAllele, AFPosteriors);
            PofFs.put(altAllele, sum);
        }
    }

    private double[] getHardyWeinbergValues(double f) {
        double[] values = hardyWeinbergValueCache.get(f);

        // if it hasn't been calculated yet, do so now
        if ( values == null ) {
            // p^2, 2pq, q^2
            values = new double[] { Math.pow(1.0 - f, 2), 2.0 * (1.0 - f) * f, Math.pow(f, 2) };
            hardyWeinbergValueCache.put(f, values);
        }

        return values;
    }

    private static double findMaxEntry(double[] array, boolean skipZeroIndex) {
        int index = skipZeroIndex ? 1 : 0;
        double max = array[index++];
        for (; index < array.length; index++) {
            if ( array[index] > max )
                max = array[index];
        }
        return max;
    }


    private enum GenotypeType { REF, HET, HOM }

    /**
     * A class for the allele-specific GL posteriors
    */
    private class AlleleSpecificGenotypeLikelihoods {
        // mapping from alt allele base to posteriors
        private HashMap<Character, double[]> alleleGLs = new HashMap<Character, double[]>();

        AlleleSpecificGenotypeLikelihoods(char ref, GenotypeLikelihoods GL) {
            double[] posteriors = GL.getPosteriors();

            // get the ref data
            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
            double refPosterior = posteriors[refGenotype.ordinal()];
            String refStr = String.valueOf(ref);

            for ( char base : BaseUtils.BASES ) {
                if ( base == ref )
                    continue;

                // get hom var and het data
                DiploidGenotype hetGenotype = ref < base ? DiploidGenotype.valueOf(refStr + String.valueOf(base)) : DiploidGenotype.valueOf(String.valueOf(base) + refStr);
                DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(base);
                double hetPosterior = posteriors[hetGenotype.ordinal()];
                double homPosterior = posteriors[homGenotype.ordinal()];

                // for precision purposes, we need to add (or really subtract, since everything is negative)
                // the largest posterior value from all entries so that numbers don't get too small
                double maxValue = Math.max(Math.max(refPosterior, hetPosterior), homPosterior);
                double normalizedRefPosterior = Math.pow(10, refPosterior - maxValue);
                double normalizedHetPosterior = Math.pow(10, hetPosterior - maxValue);
                double normalizedHomPosterior = Math.pow(10, homPosterior - maxValue);

                // normalize
                double sum = normalizedRefPosterior + normalizedHetPosterior + normalizedHomPosterior;
                normalizedRefPosterior /= sum;
                normalizedHetPosterior /= sum;
                normalizedHomPosterior /= sum;

                double[] altPosteriors = new double[] { normalizedRefPosterior, normalizedHetPosterior, normalizedHomPosterior };
                alleleGLs.put(base, altPosteriors);
            }
        }

        public double[] getNormalizedPosteriors(char altAllele) {
            return alleleGLs.get(altAllele);
        }
    }
}