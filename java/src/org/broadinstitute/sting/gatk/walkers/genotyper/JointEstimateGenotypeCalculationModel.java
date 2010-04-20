package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;

import java.util.*;

public abstract class JointEstimateGenotypeCalculationModel extends GenotypeCalculationModel {

    // for use in optimizing the P(D|AF) calculations:
    // how much off from the max likelihoods do we need to be before we can quit calculating?
    protected static final Double LOG10_OPTIMIZATION_EPSILON = 8.0;
    protected static final Double VALUE_NOT_CALCULATED = -1.0 * Double.MAX_VALUE;
    private int minAlleleFrequencyToTest;

    // because the null allele frequencies are constant for a given N,
    // we cache the results to avoid having to recompute everything
    private HashMap<Integer, double[]> nullAlleleFrequencyCache = new HashMap<Integer, double[]>();

    // because the Hardy-Weinberg values for a given frequency are constant,
    // we cache the results to avoid having to recompute everything
    // private HashMap<Double, double[]> hardyWeinbergValueCache = new HashMap<Double, double[]>();

    // the allele frequency priors
    protected double[] log10AlleleFrequencyPriors;

    // the allele frequency posteriors and P(f>0) for each alternate allele
    protected double[][] alleleFrequencyPosteriors = new double[BaseUtils.BASES.length][];
    protected double[][] log10PofDgivenAFi = new double[BaseUtils.BASES.length][];
    protected double[] PofFs = new double[BaseUtils.BASES.length];

    // the alternate allele with the largest sum of quality scores
    protected Character bestAlternateAllele = null;


    protected JointEstimateGenotypeCalculationModel() {}

    public VariantCallContext callExtendedLocus(RefMetaDataTracker tracker, char[] ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
        return null;
    }
    
    public VariantCallContext callLocus(RefMetaDataTracker tracker, char ref, GenomeLoc loc, Map<String, StratifiedAlignmentContext> contexts, DiploidGenotypePriors priors) {
        int numSamples = getNSamples(contexts);
        int frequencyEstimationPoints = (2 * numSamples) + 1;  // (add 1 for allele frequency of zero)

        // reset the optimization value
        minAlleleFrequencyToTest = 0;

        // find the alternate allele with the largest sum of quality scores
        initializeBestAlternateAllele(ref, contexts);

        // if there are no non-ref bases...
        if ( bestAlternateAllele == null ) {
            // if we don't want all bases, then we don't need to calculate genotype likelihoods
            if ( !ALL_BASE_MODE && !GENOTYPE_MODE ) {
                VariantCallContext vcc = new VariantCallContext(false);
                estimateReferenceConfidence(vcc, contexts, DiploidGenotypePriors.HUMAN_HETEROZYGOSITY, false);
                return vcc;
            }
            // otherwise, choose any alternate allele (it doesn't really matter)
            bestAlternateAllele = (ref != 'A' ? 'A' : 'C');
        }

        // calculate likelihoods if there are non-ref bases
        initializeAlleleFrequencies(frequencyEstimationPoints);

        initialize(ref, contexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE);
        calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints, contexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE);
        calculatePofFs(ref, frequencyEstimationPoints);

        // print out stats if we have a writer
        if ( verboseWriter != null )
            printAlleleFrequencyData(ref, loc, frequencyEstimationPoints);

        VariantCallContext vcc = createCalls(tracker, ref, contexts, loc, frequencyEstimationPoints);

        // technically, at this point our confidence in a reference call isn't accurately
        //  estimated because it didn't take into account samples with no data
        if ( vcc.vc == null )
            estimateReferenceConfidence(vcc, contexts, DiploidGenotypePriors.HUMAN_HETEROZYGOSITY, true);
        return vcc;
    }

    protected int getMinAlleleFrequencyToTest() {
        return minAlleleFrequencyToTest;
    }

    protected int getNSamples(Map<String, StratifiedAlignmentContext> contexts) {
        return contexts.size();
    }

    protected void initializeBestAlternateAllele(char ref, Map<String, StratifiedAlignmentContext> contexts) {
        int[] qualCounts = new int[4];

        for ( String sample : contexts.keySet() ) {
            AlignmentContext context = contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

            // calculate the sum of quality scores for each base
            ReadBackedPileup pileup = context.getBasePileup();
            for ( PileupElement p : pileup ) {
                // ignore deletions
                if ( p.isDeletion() )
                    continue;

                int index = BaseUtils.simpleBaseToBaseIndex((char)p.getBase());
                if ( index >= 0 )
                    qualCounts[index] += p.getQual();
            }
        }

        // set the non-ref base with maximum quality score sum
        int maxCount = 0;
        bestAlternateAllele = null;
        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;
            int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
            if ( qualCounts[index] > maxCount ) {
                maxCount = qualCounts[index];
                bestAlternateAllele = altAllele;
            }
        }
    }

    protected void initialize(char ref, Map<String, StratifiedAlignmentContext> contexts, StratifiedAlignmentContext.StratifiedContextType contextType) {
        // by default, no initialization is done
        return;
    }

    protected void initializeAlleleFrequencies(int frequencyEstimationPoints) {
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

    private void estimateReferenceConfidence(VariantCallContext vcc, Map<String, StratifiedAlignmentContext> contexts, double theta, boolean ignoreCoveredSamples) {

        double P_of_ref = 1.0;

        // use the AF=0 prob if it's calculated
        if ( ignoreCoveredSamples )
            P_of_ref = 1.0 - PofFs[BaseUtils.simpleBaseToBaseIndex(bestAlternateAllele)];

        // for each sample that we haven't examined yet
        for ( String sample : samples ) {
            boolean isCovered = contexts.containsKey(sample);
            if ( ignoreCoveredSamples && isCovered )
                continue;

            int depth = isCovered ? contexts.get(sample).getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().size() : 0;
            P_of_ref *= 1.0 - (theta / 2.0) * MathUtils.binomialProbability(0, depth, 0.5);
        }

        vcc.confidentlyCalled = QualityUtils.phredScaleErrorRate(1.0 - P_of_ref) >= CONFIDENCE_THRESHOLD;
    }

    protected void calculateAlleleFrequencyPosteriors(char ref, int frequencyEstimationPoints, Map<String, StratifiedAlignmentContext> contexts, StratifiedAlignmentContext.StratifiedContextType contextType) {

        // initialization
        for ( char altAllele : BaseUtils.BASES ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
            alleleFrequencyPosteriors[baseIndex] = new double[frequencyEstimationPoints];
            log10PofDgivenAFi[baseIndex] = new double[frequencyEstimationPoints];
        }

        // only calculate for the most likely alternate allele
        calculatelog10PofDgivenAFforAllF(ref, bestAlternateAllele, frequencyEstimationPoints-1, contexts, contextType);
    }

    /********************************************************************************/
    /*** One or both of the following methods should be overloaded in subclasses  ***/
    /***    so that the correct calculation is made for PofDgivenAFi              ***/
    /********************************************************************************/

    /**
     * @param ref              the ref base
     * @param alt              the alt base
     * @param numFrequencies   total number of allele frequencies (2N)
     * @param contexts         stratified alignment contexts
     * @param contextType      which stratification to use
     */
    protected void calculatelog10PofDgivenAFforAllF(char ref, char alt, int numFrequencies, Map<String, StratifiedAlignmentContext> contexts, StratifiedAlignmentContext.StratifiedContextType contextType) {
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(alt);

        // for each minor allele frequency, calculate log10PofDgivenAFi
        for (int i = 0; i <= numFrequencies; i++) {
            double f = (double)i / (double)(numFrequencies);
            log10PofDgivenAFi[baseIndex][i] += calculateLog10PofDgivenAFforF(ref, alt, f, contexts, contextType);
        }
    }

    /**
     * @param ref              the ref base
     * @param alt              the alt base
     * @param f                the allele frequency
     * @param contexts         stratified alignment contexts
     * @param contextType      which stratification to use
     *
     * @return value of PofDgivenAF for allele frequency f
     */
    protected double calculateLog10PofDgivenAFforF(char ref, char alt, double f, Map<String, StratifiedAlignmentContext> contexts, StratifiedAlignmentContext.StratifiedContextType contextType) {
        return 0.0;
    }

    /********************************************************************************/

    /**
     * @param freqI             allele frequency I
     * @param numFrequencies    total number of allele frequencies
     * @param altBaseIndex      the index of the alternate allele
     */
    protected void ignoreAlleleFrequenciesAboveI(int freqI, int numFrequencies, int altBaseIndex) {
        while ( ++freqI <= numFrequencies )
            log10PofDgivenAFi[altBaseIndex][freqI] = VALUE_NOT_CALCULATED;
    }

    /**
     * @param ref                         the ref base
     * @param frequencyEstimationPoints   number of allele frequencies
     */
    protected void calculatePofFs(char ref, int frequencyEstimationPoints) {
        // only calculate for the most likely alternate allele
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(bestAlternateAllele);

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

    /**
     * @param ref                         the ref base
     * @param loc                         the GenomeLoc
     * @param frequencyEstimationPoints   number of allele frequencies
     */
    protected void printAlleleFrequencyData(char ref, GenomeLoc loc, int frequencyEstimationPoints) {
        for (int i = 0; i < frequencyEstimationPoints; i++) {
            StringBuilder AFline = new StringBuilder("AFINFO\t");
            AFline.append(loc).append("\t");
            AFline.append(i + "/" + (frequencyEstimationPoints-1) + "\t");
            AFline.append(String.format("%.2f\t", ((float)i)/ (frequencyEstimationPoints-1)));
            AFline.append(String.format("%.8f", log10AlleleFrequencyPriors[i]) + "\t");
            for ( char altAllele : BaseUtils.BASES ) {
                if ( altAllele != ref ) {
                    int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
                    double PofDgivenAF = log10PofDgivenAFi[baseIndex][i];
                    if ( PofDgivenAF == VALUE_NOT_CALCULATED )
                        PofDgivenAF = 0.0;
                    AFline.append(String.format("%.8f\t%.8f\t", PofDgivenAF, alleleFrequencyPosteriors[baseIndex][i]));
                } else {
                    AFline.append(String.format("%.8f\t%.8f\t", -1.0, -1.0));
                }
            }
            verboseWriter.println(AFline);
        }

        for ( char altAllele : BaseUtils.BASES ) {
            if ( altAllele != ref ) {
                char base = Character.toLowerCase(altAllele);
                int baseIndex = BaseUtils.simpleBaseToBaseIndex(altAllele);
                if ( MathUtils.compareDoubles(PofFs[baseIndex], 0.0) != 0 ) {
                    double phredScaledConfidence = QualityUtils.phredScaleErrorRate(alleleFrequencyPosteriors[baseIndex][0]);
                    if ( Double.isInfinite(phredScaledConfidence) ) {
                        phredScaledConfidence = -10.0 * log10PofDgivenAFi[baseIndex][0];
                        verboseWriter.println("P(f>0)_" + base + " = 1");
                        verboseWriter.println("Qscore_" + base + " = " + phredScaledConfidence);
                        verboseWriter.println("LOD_" + base + " = " + phredScaledConfidence);
                    } else {
                        verboseWriter.println("P(f>0)_" + base + " = " + PofFs[baseIndex]);
                        verboseWriter.println("Qscore_" + base + " = " + (QualityUtils.phredScaleErrorRate(alleleFrequencyPosteriors[baseIndex][0])));
                        verboseWriter.println("LOD_" + base + " = " + (Math.log10(PofFs[baseIndex]) - Math.log10(alleleFrequencyPosteriors[baseIndex][0])));
                    }
                }
            }
        }
        verboseWriter.println();
    }

    protected Map<String, Genotype> makeGenotypeCalls(char ref, char alt, int frequency, Map<String, StratifiedAlignmentContext> contexts, GenomeLoc loc) {
        // by default, we return no genotypes
        return new HashMap<String, Genotype>();
    }    

    protected VariantCallContext createCalls(RefMetaDataTracker tracker, char ref, Map<String, StratifiedAlignmentContext> contexts, GenomeLoc loc, int frequencyEstimationPoints) {
        // only need to look at the most likely alternate allele
        int indexOfMax = BaseUtils.simpleBaseToBaseIndex(bestAlternateAllele);

        int bestAFguess = Utils.findIndexOfMaxEntry(alleleFrequencyPosteriors[indexOfMax]);
        double phredScaledConfidence;
        if ( bestAFguess != 0 ) {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(alleleFrequencyPosteriors[indexOfMax][0]);
            if ( Double.isInfinite(phredScaledConfidence) )
                phredScaledConfidence = -10.0 * log10PofDgivenAFi[indexOfMax][0];
        } else {
            phredScaledConfidence = QualityUtils.phredScaleErrorRate(PofFs[indexOfMax]);
            if ( Double.isInfinite(phredScaledConfidence) ) {
                double sum = 0.0;
                for (int i = 1; i < frequencyEstimationPoints; i++) {
                    if ( log10PofDgivenAFi[indexOfMax][i] == VALUE_NOT_CALCULATED )
                        break;
                    sum += log10PofDgivenAFi[indexOfMax][i];
                }
                phredScaledConfidence = -10.0 * sum;
            }
        }

        // return a null call if we don't pass the confidence cutoff or the most likely allele frequency is zero
        if ( !ALL_BASE_MODE && ((!GENOTYPE_MODE && bestAFguess == 0) || phredScaledConfidence < CONFIDENCE_THRESHOLD) )
            return new VariantCallContext(phredScaledConfidence >= CONFIDENCE_THRESHOLD);

        // output to beagle file if requested
        if ( beagleWriter != null ) {
            beagleWriter.print(loc);
            beagleWriter.print(' ');
            beagleWriter.print(ref);
            beagleWriter.print(' ');
            beagleWriter.print(bestAlternateAllele);
        }

        // populate the sample-specific data (output it to beagle also if requested)
        Map<String, Genotype> genotypes = makeGenotypeCalls(ref, bestAlternateAllele, bestAFguess, contexts, loc);

        // close beagle record (if requested)
        if ( beagleWriter != null )
            beagleWriter.println();

        // next, the variant context data (alleles, attributes, etc.)
        ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(new Allele(Character.toString(ref), true));
        if ( bestAFguess != 0 )
            alleles.add(new Allele(bestAlternateAllele.toString(), false));

        // *** note that calculating strand bias involves overwriting data structures, so we do that last
        HashMap<String, Object> attributes = new HashMap<String, Object>();
        if ( bestAFguess != 0 )
            attributes.put(VCFRecord.ALLELE_FREQUENCY_KEY, new Double((double)bestAFguess / (double)(frequencyEstimationPoints-1)));

        rodDbSNP dbsnp = getDbSNP(tracker);
        if ( dbsnp != null )
            attributes.put("ID", dbsnp.getRS_ID());

        if ( REPORT_SLOD ) {
            // the overall lod
            double overallLog10PofNull = log10AlleleFrequencyPriors[0] + log10PofDgivenAFi[indexOfMax][0];
            double overallLog10PofF = log10AlleleFrequencyPriors[bestAFguess] + log10PofDgivenAFi[indexOfMax][bestAFguess];
            double lod = overallLog10PofF - overallLog10PofNull;

            // set the optimization value for the subsequent strand calculations
            minAlleleFrequencyToTest = bestAFguess;

            // the forward lod
            initialize(ref, contexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD);
            calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints, contexts, StratifiedAlignmentContext.StratifiedContextType.FORWARD);
            calculatePofFs(ref, frequencyEstimationPoints);
            double forwardLog10PofNull = log10AlleleFrequencyPriors[0] + log10PofDgivenAFi[indexOfMax][0];
            double forwardLog10PofF = log10AlleleFrequencyPriors[bestAFguess] + log10PofDgivenAFi[indexOfMax][bestAFguess];

            // the reverse lod
            initialize(ref, contexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE);
            calculateAlleleFrequencyPosteriors(ref, frequencyEstimationPoints, contexts, StratifiedAlignmentContext.StratifiedContextType.REVERSE);
            calculatePofFs(ref, frequencyEstimationPoints);
            double reverseLog10PofNull = log10AlleleFrequencyPriors[0] + log10PofDgivenAFi[indexOfMax][0];
            double reverseLog10PofF = log10AlleleFrequencyPriors[bestAFguess] + log10PofDgivenAFi[indexOfMax][bestAFguess];

            double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofNull;
            double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofNull;
            //logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);

            // strand score is max bias between forward and reverse strands
            double strandScore = Math.max(forwardLod - lod, reverseLod - lod);
            // rescale by a factor of 10
            strandScore *= 10.0;
            //logger.debug(String.format("SLOD=%f", strandScore));

            attributes.put("SB", new Double(strandScore));
        }

        VariantContext vc = new VariantContext("UG_SNP_call", loc, alleles, genotypes, phredScaledConfidence/10.0, null, attributes);

        return new VariantCallContext(vc, phredScaledConfidence >= CONFIDENCE_THRESHOLD);
    }
}