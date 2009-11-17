package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.*;

public abstract class JointEstimateGenotypeCalculationModel extends GenotypeCalculationModel {

    protected JointEstimateGenotypeCalculationModel() {}

    // because the null allele frequencies are constant for a given N,
    // we cache the results to avoid having to recompute everything
    private HashMap<Integer, double[]> nullAlleleFrequencyCache = new HashMap<Integer, double[]>();

    // because the Hardy-Weinberg values for a given frequency are constant,
    // we cache the results to avoid having to recompute everything
    private HashMap<Double, double[]> hardyWeinbergValueCache = new HashMap<Double, double[]>();

    protected enum GenotypeType { REF, HET, HOM }

    // the allele frequency priors
    protected double[] log10AlleleFrequencyPriors;

    // the allele frequency posteriors and P(f>0) for each alternate allele
    protected double[][] alleleFrequencyPosteriors = new double[BaseUtils.BASES.length][];
    protected double[][] log10PofDgivenAFi = new double[BaseUtils.BASES.length][];
    protected double[] PofFs = new double[BaseUtils.BASES.length];

    
    // these need to be implemented
    protected abstract HashMap<String, AlignmentContextBySample> createContexts(AlignmentContext context);
    protected abstract void initializeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType);
    protected abstract double computeLog10PofDgivenAFi(DiploidGenotype refGenotype, DiploidGenotype hetGenotype, DiploidGenotype homGenotype, double f);
    protected abstract List<Genotype> makeGenotypeCalls(char ref, HashMap<String, AlignmentContextBySample> contexts, GenomeLoc loc);
    protected abstract int getNSamples(HashMap<String, AlignmentContextBySample> contexts);

    public Pair<List<Genotype>, GenotypeLocusData> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the context for each sample, overall and separated by strand
        HashMap<String, AlignmentContextBySample> contexts = createContexts(context);
        if ( contexts == null )
            return null;

        int numSamples = getNSamples(contexts);
        int frequencyEstimationPoints = (2 * numSamples) + 1;  // (add 1 for allele frequency of zero)

        initializeAlleleFrequencies(frequencyEstimationPoints);

        initializeLikelihoods(ref, contexts, StratifiedContext.OVERALL);
        calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints);
        calculatePofFs(ref, frequencyEstimationPoints);

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printAlleleFrequencyData(ref, context.getLocation(), frequencyEstimationPoints);

        return createCalls(tracker, ref, contexts, context.getLocation(), frequencyEstimationPoints);
   }

    private void initializeAlleleFrequencies(int frequencyEstimationPoints) {
        // set up the allele frequency priors
        log10AlleleFrequencyPriors = getNullAlleleFrequencyPriors(frequencyEstimationPoints);
    }

    protected double[] getNullAlleleFrequencyPriors(int N) {
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

    protected void calculateAlleleFrequencyPosteriors(char ref, int frequencyEstimationPoints) {

        // initialization
        for ( char altAllele : BaseUtils.BASES ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
            alleleFrequencyPosteriors[baseIndex] = new double[frequencyEstimationPoints];
            log10PofDgivenAFi[baseIndex] = new double[frequencyEstimationPoints];
        }
        DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(ref);
        String refStr = String.valueOf(ref);

        // for each alternate allele
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;

            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);

            DiploidGenotype hetGenotype = ref < altAllele ? DiploidGenotype.valueOf(refStr + String.valueOf(altAllele)) : DiploidGenotype.valueOf(String.valueOf(altAllele) + refStr);
            DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(altAllele);

            // for each minor allele frequency
            for (int i = 0; i < frequencyEstimationPoints; i++) {
                double f = (double)i / (double)(frequencyEstimationPoints-1);
                log10PofDgivenAFi[baseIndex][i] += computeLog10PofDgivenAFi(refGenotype, hetGenotype, homGenotype, f);
            }
        }
    }

    protected void calculatePofFs(char ref, int frequencyEstimationPoints) {
        // for each alternate allele
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;

            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);

            // multiply by null allele frequency priors to get AF posteriors, then normalize
            for (int i = 0; i < frequencyEstimationPoints; i++)
                alleleFrequencyPosteriors[baseIndex][i] = log10AlleleFrequencyPriors[i] + log10PofDgivenAFi[baseIndex][i];
            alleleFrequencyPosteriors[baseIndex] = MathUtils.normalizeFromLog10(alleleFrequencyPosteriors[baseIndex]);

            // calculate p(f>0)
            double sum = 0.0;
            for (int i = 1; i < frequencyEstimationPoints; i++)
                sum += alleleFrequencyPosteriors[baseIndex][i];
            PofFs[baseIndex] = Math.min(sum, 1.0); // deal with precision errors
        }
    }

    protected double[] getHardyWeinbergValues(double f) {
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

    protected void printAlleleFrequencyData(char ref, GenomeLoc loc, int frequencyEstimationPoints) {

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

    protected Pair<List<Genotype>, GenotypeLocusData> createCalls(RefMetaDataTracker tracker, char ref, HashMap<String, AlignmentContextBySample> contexts, GenomeLoc loc, int frequencyEstimationPoints) {
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
        int bestAFguess = Utils.findIndexOfMaxEntry(alleleFrequencyPosteriors[indexOfMax]);

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !ALL_BASE_MODE && (bestAFguess == 0 || phredScaledConfidence < CONFIDENCE_THRESHOLD) )
            return new Pair<List<Genotype>, GenotypeLocusData>(null, null);

        // populate the sample-specific data
        List<Genotype> calls = makeGenotypeCalls(ref, contexts, loc);

        // next, the general locus data
        // note that calculating strand bias involves overwriting data structures, so we do that last
        GenotypeLocusData locusdata = GenotypeWriterFactory.createSupportedGenotypeLocusData(OUTPUT_FORMAT, ref, loc);
        if ( locusdata != null ) {
            if ( locusdata instanceof ConfidenceBacked ) {
                ((ConfidenceBacked)locusdata).setConfidence(phredScaledConfidence);
            }
            if ( locusdata instanceof AlleleFrequencyBacked ) {
                ((AlleleFrequencyBacked)locusdata).setAlleleFrequency((double)bestAFguess / (double)(frequencyEstimationPoints-1));
                AlleleFrequencyBacked.AlleleFrequencyRange range = computeAFrange(alleleFrequencyPosteriors[indexOfMax], frequencyEstimationPoints-1, bestAFguess, ALLELE_FREQUENCY_RANGE);
                ((AlleleFrequencyBacked)locusdata).setAlleleFrequencyRange(range);
            }
            if ( locusdata instanceof IDBacked ) {
                rodDbSNP dbsnp = getDbSNP(tracker);
                if ( dbsnp != null )
                    ((IDBacked)locusdata).setID(dbsnp.getRS_ID());
            }
            if ( locusdata instanceof SLODBacked ) {
                // the overall lod
                double overallLog10PofNull = Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);
                double overallLog10PofF = Math.log10(PofFs[indexOfMax]);
                double lod = overallLog10PofF - overallLog10PofNull;

                // the forward lod
                initializeLikelihoods(ref, contexts, StratifiedContext.FORWARD);
                calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints);
                calculatePofFs(ref, frequencyEstimationPoints);
                double forwardLog10PofNull = Math.log10(alleleFrequencyPosteriors[indexOfMax][0]);
                double forwardLog10PofF = Math.log10(PofFs[indexOfMax]);

                // the reverse lod
                initializeLikelihoods(ref, contexts, StratifiedContext.REVERSE);
                calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints);
                calculatePofFs(ref, frequencyEstimationPoints);
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

    // computes the range of allele frequencies making up the given fraction of the total probability
    private static AlleleFrequencyBacked.AlleleFrequencyRange computeAFrange(double[] alleleFrequencyProbs, int N, int bestAFguess, double fraction) {
        double totalProb = alleleFrequencyProbs[bestAFguess];
        int lowIndex = bestAFguess;
        int highIndex = bestAFguess;

        // it's possible that AF=0 contains more probability than 1-fraction
        if ( alleleFrequencyProbs[0] >= (1.0 - fraction) ) {
            // in this case, the range is all possible AFs
            lowIndex = 1;
            highIndex = N;
        }
        // otherwise, find the range moving out from the best AF guess
        else {
            while ( totalProb < fraction ) {
                if ( lowIndex > 1 ) {
                    lowIndex--;
                    totalProb += alleleFrequencyProbs[lowIndex];
                }
                if ( highIndex < N ) {
                    highIndex++;
                    totalProb += alleleFrequencyProbs[highIndex];
                }
            }
        }

        return new AlleleFrequencyBacked.AlleleFrequencyRange((double)lowIndex / (double)N, (double)highIndex / (double)N, fraction);
    }
}