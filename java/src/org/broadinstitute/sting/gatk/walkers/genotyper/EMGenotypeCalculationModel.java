package org.broadinstitute.sting.gatk.walkers.genotyper;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

import java.util.*;

public abstract class EMGenotypeCalculationModel extends GenotypeCalculationModel {

    // We need to set a limit on the EM iterations in case something flukey goes on
    protected static final int MAX_EM_ITERATIONS = 6;

    // We consider the EM stable when the MAF doesn't change more than 1/10N
    protected static final double EM_STABILITY_METRIC = 0.1;

    // keep track of some metrics about our calls
    protected CallMetrics callsMetrics = new CallMetrics();


    protected EMGenotypeCalculationModel() {}

    public List<GenotypeCall> calculateGenotype(RefMetaDataTracker tracker, char ref, AlignmentContext context, DiploidGenotypePriors priors) {

        // keep track of the context for each sample, overall and separated by strand
        int[] baseCounts = new int[4];
        HashMap<String, AlignmentContextBySample> contexts = splitContextBySample(context, baseCounts);
        if ( contexts == null )
            return null;

        // debugging output
        for (int i = 0; i < 4; i++)
            logger.debug("Base count for " + BaseUtils.baseIndexToSimpleBase(i) + ": " + baseCounts[i]);

        // run the EM calculation
        EMOutput overall = runEM(ref, contexts, priors, baseCounts, StratifiedContext.OVERALL);
        double lod = overall.getPofD() - overall.getPofNull();
        logger.debug("lod=" + lod);

        // calculate strand score
        EMOutput forward = runEM(ref, contexts, priors, baseCounts, StratifiedContext.FORWARD);
        EMOutput reverse = runEM(ref, contexts, priors, baseCounts, StratifiedContext.REVERSE);
        double forwardLod = (forward.getPofD() + reverse.getPofNull()) - overall.getPofNull();
        double reverseLod = (reverse.getPofD() + forward.getPofNull()) - overall.getPofNull();
        logger.debug("forward lod=" + forwardLod + ", reverse lod=" + reverseLod);
        double strandScore = Math.max(forwardLod - lod, reverseLod - lod);

        // TODO -- finish me...

        System.out.println(String.format("LOD=%f, SLOD=%f", lod, strandScore));

        // make a call
        // List<GenotypeCall> calls
        //return new GenotypeCall(context.getLocation(), ref,gl, pileup);
        //out.addMultiSampleCall((Genotype)calls, GenotypeMetaData metadata);

        return null;
    }

    public EMOutput runEM(char ref, HashMap<String, AlignmentContextBySample> contexts, DiploidGenotypePriors priors, int[] baseCounts, StratifiedContext contextType) {

        // get initial allele frequencies
        double[] alleleFrequencies = initializeAlleleFrequencies(contexts.size(), baseCounts);
        for (int i = 0; i < alleleFrequencies.length; i++)
            logger.debug("Initial allele frequency for i=" + i + ": " + alleleFrequencies[i]);

        // get the initial genotype likelihoods
        HashMap<String, GenotypeLikelihoods> GLs = initializeGenotypeLikelihoods(ref, contexts, alleleFrequencies, priors, contextType);

        callsMetrics.nCalledBases++;
        
        // The EM loop:
        //   we want to continue until the calculation is stable, but we need some max on the number of iterations
        int iterations = 0;
        boolean EM_IS_STABLE;

        do {
            double[] newAlleleFrequencies = calculateAlleleFrequencyPosteriors(GLs);
            for (int i = 0; i < alleleFrequencies.length; i++)
                logger.debug("New allele frequency for i=" + i + ": " + newAlleleFrequencies[i]);

            applyAlleleFrequencyToGenotypeLikelihoods(GLs, newAlleleFrequencies);

            EM_IS_STABLE = isStable(alleleFrequencies, newAlleleFrequencies, contexts.size());

            alleleFrequencies = newAlleleFrequencies;

        } while ( ++iterations < MAX_EM_ITERATIONS && !EM_IS_STABLE );

        logger.debug("EM loop took " + iterations + " iterations");
        for ( String sample : GLs.keySet() )
            logger.debug("GenotypeLikelihoods for sample " + sample + ": " + GLs.get(sample).toString());

        return computePofF(ref, GLs, alleleFrequencies, contexts.size());
    }

    protected abstract double[] initializeAlleleFrequencies(int numSamplesInContext, int[] baseCounts);
    protected abstract HashMap<String, GenotypeLikelihoods> initializeGenotypeLikelihoods(char ref, HashMap<String, AlignmentContextBySample> contexts, double[] alleleFrequencies, DiploidGenotypePriors priors, StratifiedContext contextType);
    protected abstract double[] calculateAlleleFrequencyPosteriors(HashMap<String, GenotypeLikelihoods> GLs);
    protected abstract void applyAlleleFrequencyToGenotypeLikelihoods(HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies);
    protected abstract EMOutput computePofF(char ref, HashMap<String, GenotypeLikelihoods> GLs, double[] alleleFrequencies, int numSamplesInContext);


    protected boolean isStable(double[] oldAlleleFrequencies, double[] newAlleleFrequencies, int numSamplesInContext) {
        // We consider the EM stable when the MAF doesn't change more than EM_STABILITY_METRIC/N
        double EM_STABILITY_VALUE = EM_STABILITY_METRIC / (2.0 * (double)numSamplesInContext);

        // determine whether we're stable
        double AF_delta = 0.0;
        for (int i = 0; i < oldAlleleFrequencies.length; i++)
            AF_delta += Math.abs(oldAlleleFrequencies[i] - newAlleleFrequencies[i]);

        return (AF_delta < EM_STABILITY_VALUE);  
    }

    protected DiploidGenotypePriors calculateAlleleFreqBasedPriors(double[] alleleFrequencies) {
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

    /**
     * Create the mapping from sample to alignment context; also, fill in the base counts along the way.
     * Returns null iff there are too many deletions at this position.
     */
    protected HashMap<String, AlignmentContextBySample> splitContextBySample(AlignmentContext context, int[] baseCounts) {

        HashMap<String, AlignmentContextBySample> contexts = new HashMap<String, AlignmentContextBySample>();
        for (int i = 0; i < 4; i++)
            baseCounts[i] = 0;

        int deletionsInPile = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {

            // get the read and offset
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            // find the sample; special case for pools
            String sample;
            if ( POOLED_INPUT ) {
                sample = "POOL";
            } else {
                String readGroup = read.getAttribute("RG").toString(); // can't be null because those are filtered out
                sample = read.getHeader().getReadGroup(readGroup).getSample();
            }

            // create a new context object if this is the first time we're seeing a read for this sample
            AlignmentContextBySample myContext = contexts.get(sample);
            if ( myContext == null ) {
                myContext = new AlignmentContextBySample(context.getLocation());
                contexts.put(sample, myContext);
            }

            // check for deletions
            if ( offset == -1 ) {
                // are there too many deletions in the pileup?
                if ( ++deletionsInPile > maxDeletionsInPileup )
                    return null;
            } else {
                // add bad bases to the context (for DoC calculations), but don't count them
                int baseIndex = BaseUtils.simpleBaseToBaseIndex(read.getReadString().charAt(offset));
                if ( baseIndex != -1 )
                    baseCounts[baseIndex]++;
            }

            // add the read to this sample's context
            myContext.add(read, offset);
        }

        return contexts;
    }


    /**
     * A class to keep track of some basic metrics about our calls
    */
    protected class CallMetrics {
        long nConfidentCalls = 0;
        long nNonConfidentCalls = 0;
        long nCalledBases = 0;

        CallMetrics() {}

        public String toString() {
            return String.format("UG: %d confident and %d non-confident calls were made at %d bases",
                    nConfidentCalls, nNonConfidentCalls, nCalledBases);
        }
    }


    /**
     * A class to keep track of the EM output
    */
    protected class EMOutput {
        private double pD, pNull, pF;
        private HashMap<String, GenotypeLikelihoods> GLs;

        EMOutput(double pD, double pNull, double pF, HashMap<String, GenotypeLikelihoods> GLs) {
            this.pD = pD;
            this.pNull = pNull;
            this.pF = pF;
            this.GLs = GLs;
        }

        public double getPofD() { return pD; }
        public double getPofNull() { return pNull; }
        public double getPofF() { return pF; }
        public HashMap<String, GenotypeLikelihoods> getGenotypeLikelihoods() { return GLs; }
    }


    /**
     * A class to keep track of the alignment context observed for a given sample.
     * we currently store the overall context and strand-stratified sets,
     * but any arbitrary stratification could be added.
    */
    protected enum StratifiedContext { OVERALL, FORWARD, REVERSE }
    protected class AlignmentContextBySample {

        private AlignmentContext overall = null;
        private AlignmentContext forward = null;
        private AlignmentContext reverse = null;
        private GenomeLoc loc;

        private ArrayList<SAMRecord> allReads = new ArrayList<SAMRecord>();
        private ArrayList<SAMRecord> forwardReads = new ArrayList<SAMRecord>();
        private ArrayList<SAMRecord> reverseReads = new ArrayList<SAMRecord>();

        private ArrayList<Integer> allOffsets = new ArrayList<Integer>();
        private ArrayList<Integer> forwardOffsets = new ArrayList<Integer>();
        private ArrayList<Integer> reverseOffsets = new ArrayList<Integer>();


        AlignmentContextBySample(GenomeLoc loc) {
            this.loc = loc;
        }

        public AlignmentContext getContext(StratifiedContext context) {
            switch ( context ) {
                case OVERALL: return getOverallContext();
                case FORWARD: return getForwardContext();
                case REVERSE: return getReverseContext();
            }
            return null;
        }

        private AlignmentContext getOverallContext() {
            if ( overall == null )
                overall = new AlignmentContext(loc, allReads, allOffsets);
            return overall;
        }

        private AlignmentContext getForwardContext() {
            if ( forward == null )
                forward = new AlignmentContext(loc, forwardReads, forwardOffsets);
            return forward;
        }

        private AlignmentContext getReverseContext() {
            if ( reverse == null )
                reverse = new AlignmentContext(loc, reverseReads, reverseOffsets);
            return reverse;
        }

        public void add(SAMRecord read, int offset) {
            if ( read.getReadNegativeStrandFlag() ) {
                reverseReads.add(read);
                reverseOffsets.add(offset);
            } else {
                forwardReads.add(read);
                forwardOffsets.add(offset);
            }
            allReads.add(read);
            allOffsets.add(offset);
         }
    }
}