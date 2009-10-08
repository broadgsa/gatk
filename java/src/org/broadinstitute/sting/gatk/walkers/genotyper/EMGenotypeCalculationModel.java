package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.*;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

public class EMGenotypeCalculationModel extends GenotypeCalculationModel {

    // We need to set a limit on the EM iterations in case something flukey goes on
    private static final int MAX_EM_ITERATIONS = 6;

    // We consider the EM stable when the MAF doesn't change more than 1/10N
    private static final double EM_STABILITY_METRIC = 0.1;

    // keep track of some metrics about our calls
    private CallMetrics callsMetrics = new CallMetrics();


    protected EMGenotypeCalculationModel() {}

    public boolean calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the GenotypeLikelihoods for each sample, separated by strand
        HashMap<String, UnifiedGenotypeLikelihoods> GLs = new HashMap<String, UnifiedGenotypeLikelihoods>();

        // keep track of the total counts of each base in the pileup
        int[] baseCounts = new int[4];
        int totalObservedBases = 0;

        // First order of business: create the initial GenotypeLikelihoods objects.
        // Also, confirm that there aren't too many deletions in the pileup.
        int deletionsInPile = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {

            // get the read and offset
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            // skip deletions
            if ( offset == -1 ) {
                // are there too many deletions in the pileup?
                if ( ++deletionsInPile > maxDeletionsInPileup )
                    return false;
                continue;
            }

            // get the base; don't use bad bases
            char base = read.getReadString().charAt(offset);
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            if ( baseIndex == -1 )
                continue;

            // create the GL holder object if this is the first time we're seeing a base for this sample
            String readGroup = read.getAttribute("RG").toString(); // can't be null because those are filtered out
            String sample = read.getHeader().getReadGroup(readGroup).getSample();
            UnifiedGenotypeLikelihoods myGLs = GLs.get(sample);
            if ( myGLs == null ) {
                myGLs = new UnifiedGenotypeLikelihoods(baseModel, new DiploidGenotypePriors(), defaultPlatform, VERBOSE);
                GLs.put(sample, myGLs);
            }

            // assign the base to the appropriate strand
            myGLs.add(base, read, offset);

            // update the base counts
            baseCounts[baseIndex]++;
            totalObservedBases++;
        }

        // for now, we need to special-case single sample mode
        if ( samples.size() == 1 ) {
            UnifiedGenotypeLikelihoods UGL = GLs.get(samples.iterator().next());
            // if there were no good bases, the likelihoods object wouldn't exist
            if ( UGL == null )
                return false;

            callsMetrics.nCalledBases++;
            UGL.setPriors(priors);
            SSGenotypeCall call = new SSGenotypeCall(context.getLocation(), ref, UGL.getGenotypeLikelihoods(), new ReadBackedPileup(ref, context));

            if ( GENOTYPE_MODE || call.isVariant(call.getReference()) ) {
                double confidence = (GENOTYPE_MODE ? call.getNegLog10PError() : call.toVariation().getNegLog10PError());
                if ( confidence >= LOD_THRESHOLD ) {
                    callsMetrics.nConfidentCalls++;
                    out.addGenotypeCall(call);
                } else {
                    callsMetrics.nNonConfidentCalls++;
                }
            }
            return true;
        }

        callsMetrics.nCalledBases++;

        // Next, we need to create initial allele frequencies.
        // An intelligent guess would be the observed base ratios (plus some small number to account for sampling issues).

        double[] alleleFrequencies = new double[4];
        for (int i = 0; i < 4; i++)
            alleleFrequencies[i] = ((double)baseCounts[i] + DiploidGenotypePriors.HUMAN_HETEROZYGOSITY) / (double)totalObservedBases;
        DiploidGenotypePriors AFPriors = calculateAlleleFreqBasedPriors(alleleFrequencies);

        for ( UnifiedGenotypeLikelihoods UGL : GLs.values() )
            UGL.setPriors(AFPriors);

        // debugging output
        for (int i = 0; i < 4; i++)
            logger.debug("Base count and initial allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + baseCounts[i] + ", " + alleleFrequencies[i]);

        
        // The EM loop:
        //   we want to continue until the calculation is stable, but we need some max on the number of iterations
        int iterations = 0;
        // We consider the EM stable when the MAF doesn't change more than EM_STABILITY_METRIC/N
        double EM_STABILITY_VALUE = EM_STABILITY_METRIC / (2.0 * (double)GLs.size());
        boolean isStable = false;

        while ( !isStable && iterations < MAX_EM_ITERATIONS ) {

            // calculate the posterior-weighted allele frequencies and modify the priors accordingly
            double[] newAlleleFrequencies = getPosteriorWeightedFrequencies(GLs.values());
            AFPriors = calculateAlleleFreqBasedPriors(newAlleleFrequencies);
            for ( UnifiedGenotypeLikelihoods UGL : GLs.values() )
                UGL.setPriors(AFPriors);

            // determine whether we're stable
            double AF_delta = 0.0;
            for (int i = 0; i < 4; i++) {
                AF_delta += Math.abs(alleleFrequencies[i] - newAlleleFrequencies[i]);
                logger.debug("Previous allele frequency for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + alleleFrequencies[i] + ", vs. new frequency: " + newAlleleFrequencies[i]);
            }

            isStable = AF_delta < EM_STABILITY_VALUE;
            iterations++;

            alleleFrequencies = newAlleleFrequencies;
        }

        logger.debug("EM loop took " + iterations + " iterations");

        if (true)
            throw new RuntimeException("DEBUGGING");

        // compute actual priors: theta / MAF
        double pF = computeThetaOverMAF(ref, alleleFrequencies, GLs.size());

        // the posteriors from the EM loop are the population likelihoods here
        // TODO -- finish me


/*
        ClassicGenotypeLikelihoods[] G = em_result.genotype_likelihoods;
        pD = Compute_pD(G, em_result.sample_weights);
        pNull = Compute_pNull(contexts, em_result.sample_weights);

            // Apply p(f).
            double pVar = 0.0;
            for (int i = 1; i < em_result.EM_N; i++) { pVar += THETA/(double)i; }

            double p0 = Math.log10(1 - pVar);
            double pF;

            double MAF = Compute_alt_freq(ref, em_result.allele_likelihoods);

            if (MAF < 1/(2.0*em_result.EM_N)) { pF = p0; }
            else { pF = Math.log10(THETA/(2.0*em_result.EM_N * MAF)); }

            pD = pD + pF;
            pNull = pNull + p0;

        lod = pD - pNull;
        return lod;

*/
        
        //return new SSGenotypeCall(context.getLocation(), ref,gl, pileup);

        return true;
    }

    double[] getPosteriorWeightedFrequencies(Collection<UnifiedGenotypeLikelihoods> UGLs) {
        double[] newAlleleLikelihoods = new double[4];

        for ( UnifiedGenotypeLikelihoods UGL : UGLs ) {

            // calculate the posterior weighted frequencies for this sample
            double[] personalAllelePosteriors = new double[4];
            double sum = 0;
            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                double posterior = Math.pow(10, UGL.getGenotypeLikelihoods().getPosterior(g));
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base1)] += posterior;
                personalAllelePosteriors[BaseUtils.simpleBaseToBaseIndex(g.base2)] += posterior;
                sum += 2.0 * posterior;
            }

            // normalize
            for (int i = 0; i < 4; i++)
                personalAllelePosteriors[i] /= sum;
            for (int i = 0; i < 4; i++)
                newAlleleLikelihoods[i] += personalAllelePosteriors[i];
        }

        // normalize
        double sum = 0;
        for (int i = 0; i < 4; i++)
            sum += newAlleleLikelihoods[i];
        for (int i = 0; i < 4; i++)
            newAlleleLikelihoods[i] /= sum;

        return newAlleleLikelihoods;
    }

    private DiploidGenotypePriors calculateAlleleFreqBasedPriors(double[] alleleFrequencies) {
        // convert to log-space
        double[] log10Freqs = new double[4];
        for (int i = 0; i < 4; i++)
            log10Freqs[i] = Math.log10(alleleFrequencies[i]);

        double[] alleleFreqPriors = new double[10];

        // this is the Hardy-Weinberg based allele frequency (p^2, q^2, 2pq)
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            alleleFreqPriors[g.ordinal()] = log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base1)] + log10Freqs[BaseUtils.simpleBaseToBaseIndex(g.base2)];
            // add a factor of 2 for the 2pq case
            if ( g.isHet() )
                alleleFreqPriors[g.ordinal()] += Math.log10(2);
        }

        return new DiploidGenotypePriors(alleleFreqPriors);
    }

    private double computeThetaOverMAF(char ref, double[] alleleFrequencies, int samplesAtLocus) {
        double MAF;
        Integer[] sortedIndexes = Utils.SortPermutation(alleleFrequencies);
        if ( sortedIndexes[3] != BaseUtils.simpleBaseToBaseIndex(ref) )
            MAF = alleleFrequencies[sortedIndexes[3]];
        else
            MAF = alleleFrequencies[sortedIndexes[2]];

        int expectedChrs = (int)(2.0 * (double)samplesAtLocus * MAF);

        // TODO -- we need to use the priors from the UnifiedArgumentCollection
        return DiploidGenotypePriors.HUMAN_HETEROZYGOSITY / (double)expectedChrs;
    }


    /**
     * A class to keep track of some basic metrics about our calls
    */
    private class CallMetrics {
        long nConfidentCalls = 0;
        long nNonConfidentCalls = 0;
        long nCalledBases = 0;

        CallMetrics() {}

        public String toString() {
            return String.format("SSG: %d confident and %d non-confident calls were made at %d bases",
                    nConfidentCalls, nNonConfidentCalls, nCalledBases);
        }
    }
}