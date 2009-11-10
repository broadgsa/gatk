package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

public class JointEstimateGenotypeCalculationModel extends GenotypeCalculationModel {

    protected JointEstimateGenotypeCalculationModel() {}

    // because the null allele frequencies are constant for a given N,
    // we cache the results to avoid having to recompute everything
    private HashMap<Integer, double[]> nullAlleleFrequencyCache = new HashMap<Integer, double[]>();

    // because the Hardy-Weinberg values for a given frequency are constant,
    // we cache the results to avoid having to recompute everything
    private HashMap<Double, double[]> hardyWeinbergValueCache = new HashMap<Double, double[]>();

    // the allele frequency priors
    private double[] log10AlleleFrequencyPriors;

    // the allele frequency posteriors and P(f>0) for each alternate allele
    private double[][] alleleFrequencyPosteriors = new double[BaseUtils.BASES.length][];
    private double[][] log10PofDgivenAFi = new double[BaseUtils.BASES.length][];
    private double[] PofFs = new double[BaseUtils.BASES.length];

    // the minimum and actual number of points in our allele frequency estimation
    private static final int MIN_ESTIMATION_POINTS = 100;
    private int frequencyEstimationPoints;

    // the GenotypeLikelihoods map
    private HashMap<String, GenotypeLikelihoods> GLs = new HashMap<String, GenotypeLikelihoods>();

    private enum GenotypeType { REF, HET, HOM }


    public Pair<List<Genotype>, GenotypeLocusData> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the context for each sample, overall and separated by strand
        HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context);
        if ( contexts == null )
            return null;

        initializeAlleleFrequencies(contexts.size());

        // run joint estimation for the full GL contexts
        initializeGenotypeLikelihoods(ref, contexts, StratifiedContext.OVERALL);
        calculateAlleleFrequencyPosteriors(ref, context.getLocation());
        return createCalls(ref, contexts, context.getLocation());
    }

    private void initializeAlleleFrequencies(int numSamples) {

        // calculate the number of estimation points to use:
        // it's either MIN_ESTIMATION_POINTS or 2N if that's larger
        // (add 1 for allele frequency of zero)
        frequencyEstimationPoints = Math.max(MIN_ESTIMATION_POINTS, 2 * numSamples) + 1;

        // set up the allele frequency priors
        log10AlleleFrequencyPriors = getNullAlleleFrequencyPriors(frequencyEstimationPoints);
    }

    private double[] getNullAlleleFrequencyPriors(int N) {
        double[] AFs = nullAlleleFrequencyCache.get(N);

        // if it hasn't been calculated yet, do so now
        if ( AFs == null ) {

            // calculate sum(1/i)
            double sigma_1_over_I = 0.0;
            for (int i = 1; i < N; i++)
                sigma_1_over_I += 1.0 / (double)i;

            // delta = theta / sum(1/i)
            double delta = heterozygosity / sigma_1_over_I;

            // calculate the null allele frequencies for 1-N
            AFs = new double[N];
            double sum = 0.0;
            for (int i = 1; i < N; i++) {
                double value = delta / (double)i;
                AFs[i] = Math.log10(value);
                sum += value;
            }

            // null frequency for AF=0 is (1 - sum(all other frequencies))
            AFs[0] = Math.log10(1.0 - sum);

            nullAlleleFrequencyCache.put(N, AFs);
        }

        return AFs;
    }

    private void initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        GLs.clear();

        // use flat priors for GLs
        DiploidGenotypePriors priors = new DiploidGenotypePriors();

        for ( String sample : contexts.keySet() ) {
            AlignmentContextBySample context = contexts.get(sample);
            ReadBackedPileup pileup = new ReadBackedPileup(ref, context.getContext(contextType));

            // create the GenotypeLikelihoods object
            GenotypeLikelihoods GL = new GenotypeLikelihoods(baseModel, priors, defaultPlatform);
            GL.add(pileup, true);

            GLs.put(sample, GL);
        }
    }

    private void calculateAlleleFrequencyPosteriors(char ref, GenomeLoc verboseLocation) {

        // initialization
        for ( char altAllele : BaseUtils.BASES ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
            alleleFrequencyPosteriors[baseIndex] = new double[frequencyEstimationPoints];
            log10PofDgivenAFi[baseIndex] = new double[frequencyEstimationPoints];
        }
        DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
        
        // for each minor allele frequency
        for (int i = 0; i < frequencyEstimationPoints; i++) {
            double f = (double)i / (double)(frequencyEstimationPoints-1);

            // for each sample
            for ( GenotypeLikelihoods GL : GLs.values() ) {

                double[] posteriors = GL.getPosteriors();

                // get the ref data
                double refPosterior = posteriors[refGenotype.ordinal()];
                String refStr = String.valueOf(ref);

                // for each alternate allele
                for ( char altAllele : BaseUtils.BASES ) {
                    if ( altAllele == ref )
                        continue;

                    int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);

                    DiploidGenotype hetGenotype = ref < altAllele ? DiploidGenotype.valueOf(refStr + String.valueOf(altAllele)) : DiploidGenotype.valueOf(String.valueOf(altAllele) + refStr);
                    DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(altAllele);

                    double[] allelePosteriors = new double[] { refPosterior, posteriors[hetGenotype.ordinal()], posteriors[homGenotype.ordinal()] };
                    normalizeFromLog10(allelePosteriors);
                    //logger.debug("Normalized posteriors for " + altAllele + ": " + allelePosteriors[0] + " " + allelePosteriors[1] + " " + allelePosteriors[2]);

                    // calculate the posterior weighted frequencies
                    double[] HWvalues = getHardyWeinbergValues(f);
                    double PofDgivenAFi = 0.0;
                    PofDgivenAFi += HWvalues[GenotypeType.REF.ordinal()] * allelePosteriors[GenotypeType.REF.ordinal()];
                    PofDgivenAFi += HWvalues[GenotypeType.HET.ordinal()] * allelePosteriors[GenotypeType.HET.ordinal()];
                    PofDgivenAFi += HWvalues[GenotypeType.HOM.ordinal()] * allelePosteriors[GenotypeType.HOM.ordinal()];
                    log10PofDgivenAFi[baseIndex][i] += Math.log10(PofDgivenAFi);
                }
            }
        }

        // for each alternate allele
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;

            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);

            // multiply by null allele frequency priors to get AF posteriors, then normalize
            for (int i = 0; i < frequencyEstimationPoints; i++)
                alleleFrequencyPosteriors[baseIndex][i] = log10AlleleFrequencyPriors[i] + log10PofDgivenAFi[baseIndex][i];
            normalizeFromLog10(alleleFrequencyPosteriors[baseIndex]);

            // calculate p(f>0)
            double sum = 0.0;
            for (int i = 1; i < frequencyEstimationPoints; i++)
                sum += alleleFrequencyPosteriors[baseIndex][i];
            PofFs[baseIndex] = Math.min(sum, 1.0); // deal with precision errors
        }

        // print out stats if we have a position and a writer
        if ( verboseLocation != null && verboseWriter != null )
            printAlleleFrequencyData(ref, verboseLocation);
    }

    private double[] getHardyWeinbergValues(double f) {
        double[] HWvalues = hardyWeinbergValueCache.get(f);

        // if it hasn't been calculated yet, do so now
        if ( HWvalues == null ) {

            // create Hardy-Weinberg based allele frequencies (p^2, 2pq, q^2) converted to log-space
            double p = 1.0 - f;
            double q = f;

            // allele frequencies don't actually equal 0...
            if ( MathUtils.compareDoubles(q, 0.0) == 0 ) {
                q = MINIMUM_ALLELE_FREQUENCY;
                p -= MINIMUM_ALLELE_FREQUENCY;
            } else if ( MathUtils.compareDoubles(p, 0.0) == 0 ) {
                p = MINIMUM_ALLELE_FREQUENCY;
                q -= MINIMUM_ALLELE_FREQUENCY;
            }

            HWvalues = new double[] { Math.pow(p, 2), 2.0 * p * q, Math.pow(q, 2) };
            hardyWeinbergValueCache.put(f, HWvalues);
        }

        return HWvalues;
    }

    private void normalizeFromLog10(double[] array) {
        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = findMaxEntry(array).first;
        for (int i = 0; i < array.length; i++)
            array[i] = Math.pow(10, array[i] - maxValue);

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.length; i++)
            sum += array[i];
        for (int i = 0; i < array.length; i++)
            array[i] /= sum;
    }

    // returns the maximum value in the array and its index
    private static Pair<Double, Integer> findMaxEntry(double[] array) {
        int index = 0;
        double max = array[0];
        for (int i = 1; i < array.length; i++) {
            if ( array[i] > max ) {
                max = array[i];
                index = i;
            }
        }
        return new Pair<Double, Integer>(max, index);
    }

    private void printAlleleFrequencyData(char ref, GenomeLoc loc) {

        verboseWriter.println("Location=" + loc + ", ref=" + ref);
        StringBuilder header = new StringBuilder("MAF\tNullAFpriors\t");
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele != ref ) {
                char base = Character.toLowerCase(altAllele);
                header.append("P(D|AF)_" + base + "\t");
                header.append("PosteriorAF_" + base + "\t");
            }
        }
        verboseWriter.println(header);

        for (int i = 0; i < frequencyEstimationPoints; i++) {
            StringBuilder AFline = new StringBuilder(i + "/" + (frequencyEstimationPoints-1) + "\t" + String.format("%.8f", log10AlleleFrequencyPriors[i]) + "\t");
            for ( char altAllele : BaseUtils.BASES ) {
                if ( altAllele != ref ) {
                    int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
                    AFline.append(String.format("%.8f\t%.8f\t", log10PofDgivenAFi[baseIndex][i], alleleFrequencyPosteriors[baseIndex][i]));
                }
            }
            verboseWriter.println(AFline);
        }

        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele != ref ) {
                char base = Character.toLowerCase(altAllele);
                int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
                verboseWriter.println("P(f>0)_" + base + " = " + Math.log10(PofFs[baseIndex]));
                verboseWriter.println("Qscore_" + base + " = " + (-10.0 * Math.log10(alleleFrequencyPosteriors[baseIndex][0])));
                verboseWriter.println("LOD_" + base + " = " + (Math.log10(PofFs[baseIndex]) - Math.log10(alleleFrequencyPosteriors[baseIndex][0])));
            }
        }
        verboseWriter.println();
    }

    private Pair<List<Genotype>, GenotypeLocusData> createCalls(char ref, HashMap<String, AlignmentContextBySample> contexts, GenomeLoc loc) {
        // first, find the alt allele with maximum confidence
        int indexOfMax = 0;
        char baseOfMax = ref;
        double maxConfidence = Double.MIN_VALUE;
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele != ref ) {
                int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
                if ( PofFs[baseIndex] > maxConfidence ) {
                    indexOfMax = baseIndex;
                    baseOfMax = altAllele;
                    maxConfidence = PofFs[baseIndex];
                }
            }
        }

        double phredScaledConfidence = -10.0 * Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);

        // return a null call if we don't pass the confidence cutoff
        if ( !ALL_BASE_MODE && phredScaledConfidence < CONFIDENCE_THRESHOLD )
            return new Pair<List<Genotype>, GenotypeLocusData>(null, null);

        double bestAFguess = findMaxEntry(alleleFrequencyPosteriors[indexOfMax]).second / (double)(frequencyEstimationPoints-1);
        ArrayList<Genotype> calls = new ArrayList<Genotype>();

        // first, populate the sample-specific data
        for ( String sample : GLs.keySet() ) {

            // create the call
            AlignmentContext context = contexts.get(sample).getContext(StratifiedContext.OVERALL);
            Genotype call = GenotypeWriterFactory.createSupportedCall(OUTPUT_FORMAT, ref, context.getLocation());

            if ( call instanceof ReadBacked ) {
                ((ReadBacked)call).setReads(context.getReads());
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

            calls.add(call);
        }

        // next, the general locus data
        // note that calculating strand bias involves overwriting data structures, so we do that last
        GenotypeLocusData locusdata = GenotypeWriterFactory.createSupportedGenotypeLocusData(OUTPUT_FORMAT, ref, loc);
        if ( locusdata != null ) {
            if ( locusdata instanceof ConfidenceBacked ) {
                ((ConfidenceBacked)locusdata).setConfidence(phredScaledConfidence);
            }
            if ( locusdata instanceof AlleleFrequencyBacked ) {
                ((AlleleFrequencyBacked)locusdata).setAlleleFrequency(bestAFguess);
            }
            if ( locusdata instanceof AlleleBalanceBacked ) {
                ArrayList<Double> refBalances = new ArrayList<Double>();
                ArrayList<Double> onOffBalances = new ArrayList<Double>();
                ArrayList<Double> weights = new ArrayList<Double>();

                // accumulate ratios and weights
                for ( java.util.Map.Entry<String, GenotypeLikelihoods> entry : GLs.entrySet() ) {
                    // determine the best genotype
                    Integer sorted[] = Utils.SortPermutation(entry.getValue().getPosteriors());
                    DiploidGenotype bestGenotype = DiploidGenotype.values()[sorted[DiploidGenotype.values().length - 1]];

                    // we care only about het calls
                    if ( bestGenotype.isHetRef(ref) ) {

                        // make sure the alt base is our target alt
                        char altBase = bestGenotype.base1 != ref ? bestGenotype.base1 : bestGenotype.base2;
                        if ( altBase != baseOfMax )
                            continue;

                        // get the base counts at this pileup (minus deletions)
                        ReadBackedPileup pileup = new ReadBackedPileup(ref, contexts.get(entry.getKey()).getContext(StratifiedContext.OVERALL));
                        int[] counts = pileup.getBasePileupAsCounts();
                        int refCount = counts[BaseUtils.simpleBaseToBaseIndex(ref)];
                        int altCount = counts[BaseUtils.simpleBaseToBaseIndex(altBase)];
                        int totalCount = 0;
                        for (int i = 0; i < counts.length; i++)
                            totalCount += counts[i];

                        // add the entries
                        weights.add(entry.getValue().getNormalizedPosteriors()[bestGenotype.ordinal()]);
                        refBalances.add((double)refCount / (double)(refCount + altCount));
                        onOffBalances.add((double)(refCount + altCount) / (double)totalCount);
                    }
                }

                if ( weights.size() > 0 ) {
                    // normalize the weights
                    double sum = 0.0;
                    for (int i = 0; i < weights.size(); i++)
                        sum += weights.get(i);
                    for (int i = 0; i < weights.size(); i++)
                        weights.set(i, weights.get(i) / sum);

                    // calculate total weighted ratios
                    double normalizedRefRatio = 0.0;
                    double normalizedOnOffRatio = 0.0;
                    for (int i = 0; i < weights.size(); i++) {
                        normalizedRefRatio += weights.get(i) * refBalances.get(i);
                        normalizedOnOffRatio += weights.get(i) * onOffBalances.get(i);
                    }

                    ((AlleleBalanceBacked)locusdata).setRefRatio(normalizedRefRatio);
                    ((AlleleBalanceBacked)locusdata).setOnOffRatio(normalizedOnOffRatio);
                }
            }
            if ( locusdata instanceof SLODBacked ) {
                // the overall lod
                double overallLog10PofNull = Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);
                double overallLog10PofF = Math.log10(PofFs[indexOfMax]);
                double lod = overallLog10PofF - overallLog10PofNull;

                // the forward lod
                initializeGenotypeLikelihoods(ref, contexts, StratifiedContext.FORWARD);
                calculateAlleleFrequencyPosteriors(ref, loc);
                double forwardLog10PofNull = Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);
                double forwardLog10PofF = Math.log10(PofFs[indexOfMax]);

                // the reverse lod
                initializeGenotypeLikelihoods(ref, contexts, StratifiedContext.REVERSE);
                calculateAlleleFrequencyPosteriors(ref, loc);
                double reverseLog10PofNull = Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);
                double reverseLog10PofF = Math.log10(PofFs[indexOfMax]);

                double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofNull;
                double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofNull;
                //logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

                // strand score is max bias between forward and reverse strands
                double strandScore = Math.max(forwardLod - lod, reverseLod - lod);
                // rescale by a factor of 10
                strandScore *= 10.0;
                //logger.debug(String.format("SLOD=%f", strandScore));

                ((SLODBacked)locusdata).setSLOD(strandScore);
            }
        }

        return new Pair<List<Genotype>, GenotypeLocusData>(calls, locusdata);
    }
}