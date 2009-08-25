package org.broadinstitute.sting.playground.gatk.walkers.poolseq;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.ListUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import net.sf.samtools.SAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 18, 2009
 * Time: 11:57:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class CoverageAndPowerWalker extends LocusWalker<Integer, Pair<Long, Long>> {
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
        public boolean suppress_printing = false;

    @Argument(fullName="poolSize", shortName="ps", doc="Number of individuals in pool", required=true)
        public int num_individuals = 0;

    @Argument(fullName="bootStrap", shortName="bs", doc="Use a bootstrap method", required=false)
        public boolean use_bootstrap = false;

    @Argument(fullName="lodThreshold", shortName="lt", doc="Threshold for LOD score for calls")
        public double threshold = 3.0;


     
    private static final int BOOTSTRAP_ITERATIONS = 300;


    @Override
    public void initialize()
    {

        if(num_individuals <= 0)
            throw new IllegalArgumentException("Positive nonzero parameter expected for poolSize");
    }


    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
        if ( !suppress_printing )
        {
            Pair<Double,Byte> powpair = boostrapSamplingPowerCalc(context);
            out.printf("%s: %d %d %f%n", context.getLocation(), context.getReads().size(),powpair.second,powpair.first);
        }

        return context.getReads().size();
    }


    public boolean isReduceByInterval() {
        return true;
    }


    public Pair<Long, Long> reduceInit() {

        return new Pair<Long,Long>(0l,0l); }


    public Pair<Long, Long> reduce(Integer value, Pair<Long, Long> sum)
    {
        long left = value.longValue() + sum.getFirst();
        long right = sum.getSecond() + 1l;
        return new Pair<Long,Long>(left, right);
    }


    public void onTraversalDone(Pair<Long, Long> result) {
        out.printf("Average depth of coverage is: %.2f in %d total coverage over %d sites\n",
                ((double)result.getFirst() / (double)result.getSecond()), result.getFirst(), result.getSecond());
    }

    /*
     *
     * Helper Methods Below
     *
     */

    public Pair<Double,Byte> powerTheoretical(int depth, List<SAMRecord> reads, List<Integer> offsets, double snp_prop)
    {
        // get the median Q score for these reads

        byte medQ = ListUtils.getQScoreMedian(reads,offsets);
        // System.out.println("Median Q: " + medQ); //TODO: remove this line
        double p_error = QualityUtils.qualToErrorProb(medQ);

        // variable names from here on out come from the mathematics of a log likelihood test
        // with binomial probabilities, of observing k bases consistent with the same SNP
        // given that SNP occurs in one of the alleles, versus it does not occur.
        // Numerator and denominator will each have a term raised to the kth power
        // and a term raised to the (depth - kth) (or d-kth) power. Thus the names.
        // Taking a log then brings those powers out as coefficients, where they can be solved for.


        double kterm_num = Math.log10(snp_prop * (1 - p_error) + (1 - snp_prop) * (p_error/3));
        double kterm_denom = Math.log10(p_error/3);
        double dkterm_num = Math.log10(snp_prop*(p_error/3) + (1 - snp_prop) * (1 - p_error));
        double dkterm_denom = Math.log10(1 - p_error);

        int kaccept = (int) Math.ceil((threshold-((double)depth)*(dkterm_num-dkterm_denom))/(kterm_num-kterm_denom-dkterm_num+dkterm_denom));

        // we will reject the null hypothesis if we see kaccept or more SNPs, the power is the probability that this occurs
        // we can optimize this by checking to see which sum is smaller

        double pow = 0;

        if(depth - kaccept < kaccept) {// kaccept > depth/2 - calculate power as P[hits between k and depth]

            for(int k = kaccept; k < depth; k++) {

                pow += MathUtils.binomialProbabilityLog(k, depth, snp_prop);
            }

            return new Pair(pow,medQ);
        } else { // kaccept < depth/2 - calculate power as 1-P[hits between 0 and k]
            for(int k = 0; k < kaccept; k++) {
                pow += MathUtils.binomialProbabilityLog(k,depth,snp_prop);
            }

            return new Pair(1.0-pow,medQ);
        }
    }


    public Pair<Double,Byte> boostrapSamplingPowerCalc(AlignmentContext context)
    {
        /* it is assumed that there are 2*num_individuals chromosomes in the pool
         * we use a bootstrap method to calculate our power to detect a SNP with this
         * depth of coverage
         * Assume that only one of the individual chromosomes has a variant at this locus
         *
         */

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        final int depth = reads.size();
        Pair<Double,Byte> result;
        final double single_snip_proportion = 1/(2.0*num_individuals);

        if (depth <= 0) {
            result = new Pair(-1,-1);
        } else if (!use_bootstrap) { // object data from command line

            result = powerTheoretical(depth, reads, offsets, single_snip_proportion);
        } else {

            //
            // otherwise, bootstrapping occurs below
            //

            int hypothesis_rejections=0;

            for(int boot = 0; boot < BOOTSTRAP_ITERATIONS; boot++)
            {
                List<Byte> qscore_SNPs = new LinkedList<Byte>();
                List<Byte> qscore_refs = new LinkedList<Byte>();

                for (int strap = 0; strap < depth; strap++) {
                    int readOffset = randomlySelectRead(depth);
                    byte baseQuality = reads.get(readOffset).getBaseQualities()[offsets.get(readOffset)];

                    if (Math.random() < single_snip_proportion) {// occurs with probability "single_snip_proportion"
                        qscore_SNPs.add(baseQuality);
                    }
                    else { // simulating a reference read
                        qscore_refs.add(baseQuality);
                    }
                }

                if (calculateLODByQLists(qscore_SNPs,qscore_refs) >= threshold) {
                    hypothesis_rejections++;
                }

            }

            result = new Pair(((double)hypothesis_rejections)/BOOTSTRAP_ITERATIONS, ListUtils.getQScoreMedian(reads,offsets));
        }

        return result;

    }


    public int randomlySelectRead(int depth)
    {
        double qscore_selector = Math.random();
        int readspositionrandom;
        for(readspositionrandom = 1; readspositionrandom < ((double)depth * qscore_selector); readspositionrandom ++) {
            if(readspositionrandom > depth + 1) {
                throw new RuntimeException("qscore iterator exceeding possible thresholds");
            }
        }

        return readspositionrandom - 1;
    }


    public double calculateLODByQLists(List<Byte> q_snp, List<Byte> q_ref)
    {
        /* LOD score calculation
         * let Qs be the list q_snp and Qr q_ref and f be a function that
         * converts Q-scores to probability of error, then the probability ratio
         * is given by
         *
         * Product_{i=1}^{|Qs|}((1-f(Qs[i])/2n+(2n-1)f(Qs[i])/2n)*Product_{j=1}^{|Qr|}(((2n-1)(1-f(Qr[j]))+f(Qr[j]))/2n)
         * ------------------------------------------------------------------------------------------------------------
         *                        Product_{i=1}^{|Qs|}f(Qs[i])*Product_{j=1}^{|Qr|}(1-f(Qr[j]))
         *
         * since depth = |Qr|+|Qs| the binomial coefficients cancel so are not present in this formula
         *
         *
         *
         *
         * the "LOD" used here is really a log-likelihood, just the logarithm of the above
         *
         * Helper functions compute all but the final sum (log of mult/div becomes add/subtract)
         *
         */

        Pair<Double,Double> logsumSNP = qListToSumLogProbabilities(true,q_snp);
        Pair<Double,Double> logsumRef = qListToSumLogProbabilities(false,q_ref);

        return 0 - logsumSNP.first - logsumRef.first + logsumSNP.second + logsumRef.second;
    }


    public Pair<Double,Double> qListToSumLogProbabilities(boolean listRepresentsSNPObservations, List<Byte> qList)
    {
        double logProbObserveXAndSNPTrue = 0; // note "error" for SNP is observing a ref
        double logProbObserveXAndRefTrue = 0;// and "error" for ref is observing a SNP
        final double denom = 2*((double)num_individuals);

        for (byte qual : qList) {
            double p_err = QualityUtils.qualToErrorProb(qual);
            if (listRepresentsSNPObservations) {
                logProbObserveXAndSNPTrue += Math.log10((1 - p_err) / denom +((denom - 1)*p_err) / denom);
                logProbObserveXAndRefTrue += Math.log10(p_err);
            } else {
                logProbObserveXAndSNPTrue += Math.log10((denom - 1) * (1 - p_err)/denom + p_err/denom);
                logProbObserveXAndRefTrue+= Math.log10(1 -p_err);
            }
        }

        return new Pair<Double,Double>(logProbObserveXAndSNPTrue,logProbObserveXAndRefTrue);
    }



}