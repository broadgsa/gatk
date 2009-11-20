package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.util.*;

import net.sf.samtools.SAMRecord;

public class PooledCalculationModel extends JointEstimateGenotypeCalculationModel {
    private static final String POOL_SAMPLE_NAME = "POOL";

    private FourBaseProbabilities fourBaseLikelihoods = null;

    private double[] log10F = null;
    private double[] log10OneMinusF = null;

    /**
     *
     */
    protected PooledCalculationModel() {}

    /**
     *
     * @param ref
     * @param contexts
     * @param contextType
     */
    protected void initialize(char ref, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        super.initialize(ref, contexts, contextType);

        // prepare the four base vector calculator
        fourBaseLikelihoods = FourBaseProbabilitiesFactory.makeFourBaseLikelihoods(baseModel, defaultPlatform);

        // todo - can you add static initialize()
        // prepare cached values of the log10 of f and (1-f)
        if ( log10F == null ) {
            int nChromosomes = 2 * getNSamples(contexts);
            log10F = new double[nChromosomes+1];
            log10OneMinusF = new double[nChromosomes+1];
            for ( int i = 0; i < (nChromosomes+1); i++ ) {
                double f = (1.0 * i) / nChromosomes;
                log10F[i] = Math.log10(f);
                log10OneMinusF[i] = Math.log10(1-f);
            }
        }   
    }

    protected int getNSamples(HashMap<String, AlignmentContextBySample> contexts) {
        return POOL_SIZE;
    }
    
    protected HashMap<String, AlignmentContextBySample> createContexts(AlignmentContext context) {
        // for testing purposes, we may want to throw multi-samples at pooled mode,
        // so we can't use the standard splitContextBySample() method here
        AlignmentContextBySample pooledContext = new AlignmentContextBySample(context.getLocation());

        int deletionsInPileup = 0;
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int i = 0; i < reads.size(); i++) {
            // check for deletions
            int offset = offsets.get(i);
            if ( offset == -1 ) {
                // are there too many deletions in the pileup?
                if ( ++deletionsInPileup > maxDeletionsInPileup && maxDeletionsInPileup >= 0 )
                    return null;
            }

            // add the read to this sample's context
            // note that bad bases are added to the context (for DoC calculations later)
            pooledContext.add(reads.get(i), offset);
        }

        HashMap<String, AlignmentContextBySample> contexts = new HashMap<String, AlignmentContextBySample>();
        contexts.put(POOL_SAMPLE_NAME, pooledContext);
        return contexts;
    }

    //
    // todo - create cache indexed by chromosome, tech, qual, observed base, and potential polymorphism (A/C) for all
    // todo -- 16 pairs.  We can change the calculation at a locus to: for each read in pileup, calculate offset into
    // todo -- this table given observed base, read, and qual.  Then, for each i -> 0, 2n and polymorphism, run over
    // todo -- the table, summing values
    //
    protected double computeLog10PofDgivenAFi(char refArg, char altArg, double f, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        AlignmentContextBySample context = contexts.get(POOL_SAMPLE_NAME);
        ReadBackedPileup pileup = new ReadBackedPileup(refArg, context.getContext(contextType));

        double log10L = 0.0;

        int refIndex = BaseUtils.simpleBaseToBaseIndex(refArg);
        int altIndex = BaseUtils.simpleBaseToBaseIndex(altArg);

        int nChromosomes = 2 * getNSamples(contexts);
        int nAltAlleles = (int)(f * nChromosomes);
        int nRefAlleles = nChromosomes - nAltAlleles;

        double log10POfRef = log10OneMinusF[nAltAlleles];
        double log10POfAlt = log10F[nAltAlleles];
        double POfRef = Math.pow(10,log10POfRef);
        double POfAlt = Math.pow(10,log10POfAlt);

        for (int i = 0; i < pileup.getReads().size(); i++) {
            int offset = pileup.getOffsets().get(i);

            // ignore deletions
            if ( offset == -1 )
                continue;

            SAMRecord read = pileup.getReads().get(i);
            char base = (char)read.getReadBases()[offset];
            int bIndex = BaseUtils.simpleBaseToBaseIndex(base);
            byte qual = read.getBaseQualities()[offset];

            if ( qual > 0 && bIndex != -1 ) {
                FourBaseProbabilities fbl = fourBaseLikelihoods.computeLog10Likelihoods(base, qual, read, offset);
                double POfBGivenRef = fbl.getLikelihood(refIndex);
                double POfBGivenAlt = fbl.getLikelihood(altIndex);
                double P = POfRef * POfBGivenRef + POfAlt * POfBGivenAlt;
                log10L += Math.log10(P);
            }
        }

        //if ( verboseWriter != null )
        //    verboseWriter.printf("POOL_DEBUG %s %c %c %d %d %d %.2f %.2f %.2f %f%n",
        //            context.getContext(StratifiedContext.OVERALL).getLocation(),
        //            refArg, altArg, nChromosomes, nAltAlleles, nRefAlleles, f, log10POfRef, log10POfAlt, log10L);

        return log10L;
    }

/*    protected double computeLog10PofDgivenAFi_V2(char refArg, char altArg, double f, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        AlignmentContextBySample context = contexts.get(POOL_SAMPLE_NAME);
        ReadBackedPileup pileup = new ReadBackedPileup(refArg, context.getContext(contextType));

        double log10L = 0.0;

        int refIndex = BaseUtils.simpleBaseToBaseIndex(refArg);
        int altIndex = BaseUtils.simpleBaseToBaseIndex(altArg);

        int nChromosomes = 2 * getNSamples(contexts);
        int nAltAlleles = (int)(f * nChromosomes);
        int nRefAlleles = nChromosomes - nAltAlleles;

        double log10POfRef = log10OneMinusF[nAltAlleles];
        double log10POfAlt = log10F[nAltAlleles];

        for (int i = 0; i < pileup.getReads().size(); i++) {
            int offset = pileup.getOffsets().get(i);

            // ignore deletions
            if ( offset == -1 )
                continue;

            SAMRecord read = pileup.getReads().get(i);
            char base = (char)read.getReadBases()[offset];
            int bIndex = BaseUtils.simpleBaseToBaseIndex(base);
            byte qual = read.getBaseQualities()[offset];

            double log10POfNotB = Math.log10(QualityUtils.qualToErrorProb(qual));
            if ( qual > 0 && bIndex != -1 ) {
                FourBaseProbabilities fbl = fourBaseLikelihoods.computeLog10Likelihoods(base, qual, read, offset);

                double log10POfB = fbl.getLog10Likelihood(base);

                if ( bIndex == refIndex && nRefAlleles > 0 ) {
                    log10L += log10POfRef + log10POfB;
                } else if ( bIndex == altIndex && nAltAlleles > 0) {
                    log10L += log10POfAlt + log10POfB;
                } else {
                    //log10L += Math.min(log10POfRef + fbl.getLog10Likelihood(refIndex), log10POfAlt + fbl.getLog10Likelihood(altIndex));
                    log10L += log10POfNotB;
                }
            }
        }

        if ( verboseWriter != null )
            verboseWriter.printf("POOL_DEBUG %s %c %c %d %d %d %.2f %.2f %.2f %f%n",
                    context.getContext(StratifiedContext.OVERALL).getLocation(),
                    refArg, altArg, nChromosomes, nAltAlleles, nRefAlleles, f, log10POfRef, log10POfAlt, log10L);

        return log10L;
    } */

/*    protected double computeLog10PofDgivenAFi_V1(char refArg, char altArg, double f, HashMap<String, AlignmentContextBySample> contexts, StratifiedContext contextType) {
        AlignmentContextBySample context = contexts.get(POOL_SAMPLE_NAME);
        ReadBackedPileup pileup = new ReadBackedPileup(refArg, context.getContext(contextType));

        double log10L = 0.0;

        int refIndex = BaseUtils.simpleBaseToBaseIndex(refArg);
        int altIndex = BaseUtils.simpleBaseToBaseIndex(altArg);

        int nChromosomes = 2 * getNSamples(contexts);
        int nAltAlleles = (int)(f * nChromosomes);
        int nRefAlleles = nChromosomes - nAltAlleles;

        double log10POfRef = Math.log10(1 - f);
        double log10POfAlt = Math.log10(f);
        //double log10ChromChooseRef = Math.log10(Arithmetic.binomial(nChromosomes, nRefAlleles));
        //double log10ChromChooseAlt = Math.log10(Arithmetic.binomial(nChromosomes, nAltAlleles));

        byte[] bases = pileup.getBases();
        byte[] quals = pileup.getQuals();

        for ( int i = 0; i < bases.length; i++ ) {
            int bIndex = BaseUtils.simpleBaseToBaseIndex((char)bases[i]);
            byte qual = quals[i];
            if ( qual > 0 && bIndex != -1 ) {
                double log10POfB = Math.log10(QualityUtils.qualToProb(qual));
                double log10POfNotB = Math.log10(QualityUtils.qualToErrorProb(qual));
                //System.out.printf("%f %f %f %d%n", log10L, log10POfB, log10POfNotB, qual);

                if ( bIndex == refIndex && nRefAlleles > 0 ) {
                    log10L += log10POfRef + log10POfB;
                } else if ( bIndex == altIndex && nAltAlleles > 0) {
                    log10L += log10POfAlt + log10POfB;
                } else {
                    log10L += log10POfNotB;
                }
            }
        }

        if ( verboseWriter != null )
            verboseWriter.printf("POOL_DEBUG %s %c %c %d %d %d %.2f %.2f %.2f %f%n",
                    context.getContext(StratifiedContext.OVERALL).getLocation(),
                    refArg, altArg, nChromosomes, nAltAlleles, nRefAlleles, f, log10POfRef, log10POfAlt, log10L);

        return log10L;
    }*/
}