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
import org.broadinstitute.sting.playground.utils.PoolUtils;
import net.sf.samtools.SAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 18, 2009
 * Time: 11:57:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class CoverageAndPowerWalker extends LocusWalker<Pair<Integer, Integer>, Pair<Long, Long>> {
    @Argument(fullName="suppressLocusPrinting",doc="Suppress printing",required=false)
        public boolean suppress_printing = false;

    @Argument(fullName="poolSize", shortName="ps", doc="Number of individuals in pool", required=true)
        public int num_individuals = 0;

    @Argument(fullName="bootStrap", shortName="bs", doc="Use a bootstrap method", required=false)
        public boolean useBootstrap = false;

    @Argument(fullName="lodThreshold", shortName="lt", doc="Threshold for LOD score for calls")
        public double threshold = 3.0;
     
    private static final int BOOTSTRAP_ITERATIONS = 300;

    public void initialize() {

        if(num_individuals <= 0)
            throw new IllegalArgumentException("Positive nonzero parameter expected for poolSize");
    }

    public Pair<Long,Long> reduceInit() {
        if ( ! suppress_printing ) { // print header
            out.printf("%s%n",createHeaderString());
        }
        return new Pair(0l,0l);
    }

    public Pair<Long,Long> reduce(Pair<Integer,Integer> newCvg, Pair<Long,Long> prevCvg) {
        return new Pair(prevCvg.first + newCvg.first, prevCvg.second + newCvg.second);
    }

    public Pair<Integer,Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> readsByDirection = PoolUtils.splitReadsByReadDirection(context.getReads(),context.getOffsets());
        if ( ! suppress_printing) {
            Pair<double[],byte[]> powers = calculatePower(readsByDirection, useBootstrap, context);
            out.printf("%s: %d %d %d %d %d %d %f %f %f%n", context.getLocation(), readsByDirection.getFirst().getFirst().size(), readsByDirection.getFirst().getSecond().size(),
            context.getReads().size(), powers.getSecond()[0], powers.getSecond()[1], powers.getSecond()[2],
            powers.getFirst()[0], powers.getFirst()[1], powers.getFirst()[2]);
        }
        return new Pair(readsByDirection.getFirst().getFirst().size(),readsByDirection.getFirst().getSecond().size());
    }

    // helper methods

    public Pair<double[],byte[]> calculatePower(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> dirReads, boolean bootStrap, AlignmentContext context) {
        Pair<double[],byte[]> powerAndMedianQ;
        if( bootStrap ) {
            powerAndMedianQ = calculateDirectionalPowersByBootstrap(dirReads,threshold,context);
        } else {
            powerAndMedianQ = calculateDirectionalPowersTheoretical(dirReads,threshold,context);
        }

        return powerAndMedianQ;
    }

    public Pair<double[],byte[]> calculateDirectionalPowersTheoretical(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> dirReads, double thresh, AlignmentContext context) {
        byte[] medQs = getMedianQScores(dirReads, context);
        double powerForward = calculatePowerTheoretical(dirReads.getFirst().getFirst(),dirReads.getSecond().getFirst(),medQs[0],thresh);
        double powerReverse = calculatePowerTheoretical(dirReads.getFirst().getSecond(), dirReads.getSecond().getSecond(),medQs[1],thresh);
        double powerCombined = calculatePowerTheoretical(context.getReads(),context.getOffsets(),medQs[2],thresh);
        double [] powers = {powerForward,powerReverse,powerCombined};
        return new Pair(powers, medQs);
    }

    public byte[] getMedianQScores(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> dirReads, AlignmentContext context) {
        byte medQForward = getMedianQualityScore(dirReads.getFirst().getFirst(),dirReads.getSecond().getFirst());
        byte medQReverse = getMedianQualityScore(dirReads.getFirst().getSecond(), dirReads.getSecond().getSecond());
        byte medQCombined = getMedianQualityScore(context.getReads(),context.getOffsets());
        byte [] medQs = {medQForward, medQReverse, medQCombined};
        return medQs;
    }

    public byte getMedianQualityScore(List<SAMRecord> reads, List<Integer> offsets) {
        return ListUtils.getQScoreMedian(reads,offsets);
    }

    public double calculatePowerTheoretical(List<SAMRecord> reads, List<Integer> offsets, byte medQ, double thresh) {
        double p_error = QualityUtils.qualToErrorProb(medQ);
        int depth = reads.size();
        if(depth <= 0) {
            return 0;
        }
        double snpProp = this.getSNPProportion(1);

        // variable names from here on out come from the mathematics of a log likelihood test
        // with binomial probabilities, of observing k bases consistent with the same SNP
        // given that SNP occurs in one of the alleles, versus it does not occur.
        // Numerator and denominator will each have a term raised to the kth power
        // and a term raised to the (depth - kth) (or d-kth) power. Thus the names.
        // Taking a log then brings those powers out as coefficients, where they can be solved for.


        double kterm_num = Math.log10( snpProp * (1 - p_error) + (1 - snpProp) * (p_error/3) );
        double kterm_denom = Math.log10( p_error/3 );
        double dkterm_num = Math.log10( snpProp * (p_error/3) + (1 - snpProp) * (1 - p_error) );
        double dkterm_denom = Math.log10( 1 - p_error);

        int kaccept = (int) Math.ceil( ( thresh - ( (double) depth ) * ( dkterm_num - dkterm_denom ) )/( kterm_num - kterm_denom- dkterm_num + dkterm_denom ) );

        // we will reject the null hypothesis if we see kaccept or more SNPs, the power is the probability that this occurs
        // we can optimize this by checking to see which sum is smaller

        double pow = 0;

        if ( depth - kaccept < kaccept ) {// kaccept > depth/2 - calculate power as P[hits between kaccept and depth]
            return MathUtils.cumBinomialProbLog(kaccept, depth, depth, snpProp);
        } else { // kaccept < depth/2 - calculate power as 1-P[hits between 0 and kaccept]
            return 1-MathUtils.cumBinomialProbLog(0,kaccept,depth,snpProp);
        }
    }

    public double getSNPProportion(int numSNPAlleles) {
        return ((double) numSNPAlleles / (2.0 * num_individuals));
    }

    public Pair<double[],byte[]> calculateDirectionalPowersByBootstrap(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> dirReads, double thresh, AlignmentContext context) {
        double powerForward = bootstrapPower(dirReads.getFirst().getFirst(),dirReads.getSecond().getFirst(),thresh);
        double powerReverse = bootstrapPower(dirReads.getFirst().getSecond(),dirReads.getSecond().getSecond(),thresh);
        double powerCombined = bootstrapPower(context.getReads(),context.getOffsets(),thresh);
        double[] powers = {powerForward, powerReverse, powerCombined};
        return new Pair(powers,getMedianQScores(dirReads,context));
    }

    public double bootstrapPower(List<SAMRecord> reads, List<Integer> offsets, double thresh) {
        double power = 0.0;
        for ( int boot = 0; boot < BOOTSTRAP_ITERATIONS; boot++) {
            Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> snpReadsAndRefReads = coinTossPartition(reads,offsets,this.getSNPProportion(1));
            if( calculateLogLikelihoodOfSample(snpReadsAndRefReads, num_individuals) > thresh) {
                power += 1.0/BOOTSTRAP_ITERATIONS;
            }
        }

        return power;
    }

    public String createHeaderString() {
        return "Chrom:Pos  CvgF  CvgR  CvgC  QmedF  QmedR  QmedC  PowF  PowR  PowC";
    }

    // class methods

    public static Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> coinTossPartition(List<SAMRecord> reads, List<Integer> offsets, double snpProb) {
        Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> partitionedReads;
        if( reads == null || offsets == null ) {
            partitionedReads = null;
        } else {
            ArrayList<SAMRecord> snpReads = new ArrayList();
            ArrayList<SAMRecord> refReads = new ArrayList();
            ArrayList<Integer> snpOff = new ArrayList();
            ArrayList<Integer> refOff = new ArrayList();
            for(int readNo = 0; readNo < reads.size(); readNo++) {
                if ( Math.random() < snpProb) {
                    snpReads.add(reads.get(readNo));
                    snpOff.add(offsets.get(readNo));
                } else {
                    refReads.add(reads.get(readNo));
                    refOff.add(offsets.get(readNo));
                }
            }
            partitionedReads= new Pair(new Pair(snpReads,refReads), new Pair(snpOff, refOff));
        }

        return partitionedReads;
    }

    public static double calculateLogLikelihoodOfSample(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> snpReadsRefReads, int nIndivids) {
        List<Byte> qListSnps = getQList(snpReadsRefReads.getFirst().getFirst(),snpReadsRefReads.getSecond().getFirst());
        List<Byte> qListRefs = getQList(snpReadsRefReads.getFirst().getSecond(),snpReadsRefReads.getSecond().getSecond());
        Pair<Double,Double> logsumSNP = qListToSumLogProbabilities(true,qListSnps, 2.0*nIndivids);
        Pair<Double,Double> logsumRef = qListToSumLogProbabilities(false,qListRefs, 2.0*nIndivids);

        return 0.0 - logsumSNP.first - logsumRef.first + logsumSNP.second + logsumRef.second;
    }

    public static List<Byte> getQList(List<SAMRecord> reads, List<Integer> offsets) {
        List<Byte> qscores = new LinkedList();
        for(int readNo = 0; readNo < reads.size(); readNo++) {
            qscores.add(reads.get(readNo).getBaseQualities()[offsets.get(readNo)]);
        }

        return qscores;
    }

    public static Pair<Double,Double> qListToSumLogProbabilities(boolean listRepresentsSNPObservations, List<Byte> qList, double denom)
    {
        double logProbObserveXAndSNPTrue = 0; // note "error" for SNP is observing a ref
        double logProbObserveXAndRefTrue = 0;// and "error" for ref is observing a SNP

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