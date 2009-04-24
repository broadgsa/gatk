package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.MathUtils;

import net.sf.samtools.SAMRecord;

import java.util.*;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, Integer>  {
    @Argument(fullName="metrics",required=true)
    public String metricsFileName;

    @Argument(fullName="lodThreshold",shortName="lod",required=false,defaultValue="5.0")
    public Double lodThreshold;

    @Argument(fullName="fourBaseMode",required=false,defaultValue="false")
    public Boolean fourBaseMode;

    @Argument(fullName="retest",required=false,defaultValue="false")
    public Boolean retest;

    @Argument(fullName="printMetrics",required=false,defaultValue="true")
    public Boolean printMetrics;

    public AlleleMetrics metrics;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) { return true; }
    public boolean requiresReads() { return true; }
    public void initialize() {
        metrics = new AlleleMetrics(metricsFileName, lodThreshold);
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        String rodString = getRodString(tracker);

        AlleleFrequencyEstimate freq = null;
        if (fourBaseMode) 
        {
            if (retest) {
                // Compute single quality score genotype likelihoods
                freq = getOneProbAlleleFrequency(ref, context, rodString);

                if (MathUtils.compareDoubles(freq.qhat, 0.5) == 0) {
                    // If we found a putative het, recompute four-base prob genotype likelihoods
                    freq = getFourProbAlleleFrequency(ref, context, rodString, tracker);
                }
            } else {
                freq = getFourProbAlleleFrequency(ref, context, rodString, tracker);
            }
        }
        else 
        {
            // Compute single quality score genotype likelihoods
            freq = getOneProbAlleleFrequency(ref, context, rodString);
        }

        if (printMetrics) {
            if (freq != null) { metrics.nextPosition(freq, tracker); }
            metrics.printMetricsAtLocusIntervals(1000);
        }

        return freq;
    }

    private AlleleFrequencyEstimate getFourProbAlleleFrequency(char ref, LocusContext context, String rodString, RefMetaDataTracker tracker) {
        /*
         P(D)*P(G,q|D) = P(D|G,q)*P(G,q)
                       = P(D|q)*P(q|G)*P(G)

                  P(G) = { 0.999 (hom-ref), 1e-3 (het), 1e-5 (hom-nonref) }
                  
                            n
          P(D|q)P(q|G) = product [ P_i(A)*B(i, n, G) ]
                           i=1
         */

        double[][] probs = ReadBackedPileup.probDistPileup(context.getReads(), context.getOffsets());
        int refBaseIndex = BaseUtils.simpleBaseToBaseIndex(ref);
        int altBaseIndex = getMostFrequentNonRefBase(getFractionalCounts(probs), refBaseIndex);

        if (refBaseIndex >= 0 && altBaseIndex >= 0) {
            double[] obsWeights = getObservationWeights(probs, refBaseIndex, altBaseIndex);

            double[] genotypePriors   = { 0.999, 1e-3, 1e-5 };
            double[] genotypeBalances = { 0.02, 0.5,  0.98 };
            double[] reportingBalances = { 0.0, 0.5, 1.0 }; // the metrics class doesn't like it if you don't tell it these exact values

            double[] posteriors = new double[3];
            double qhat = 0.0, qstar = 0.0, lodVsRef = 0.0, pBest = Double.MIN_VALUE;

            for (int hypothesis = 0; hypothesis < 3; hypothesis++) {
                posteriors[hypothesis] = 0.0;

                for (int weightIndex = 0; weightIndex < obsWeights.length; weightIndex++) {
                    double binomialWeight = MathUtils.binomialProbability(weightIndex, probs.length, genotypeBalances[hypothesis]);
                    
                    //if (hypothesis == 0 && weightIndex < obsWeights.length/4) {
                    //    binomialWeight = 1.0;
                    //}

                    posteriors[hypothesis] += obsWeights[weightIndex]*binomialWeight;
                }
                posteriors[hypothesis] *= genotypePriors[hypothesis];

                double refLikelihood = Double.isInfinite(Math.log10(posteriors[0])) ? -10000.0 : Math.log10(posteriors[0]);

                if (posteriors[hypothesis] > pBest) {
                    qhat = reportingBalances[hypothesis];
                    qstar = reportingBalances[hypothesis];
                    lodVsRef = Math.log10(posteriors[hypothesis]) - refLikelihood;
                    pBest = posteriors[hypothesis];
                }
            }

            double pRef = posteriors[0];

            double[] sposteriors = Arrays.copyOf(posteriors, posteriors.length);
            Arrays.sort(sposteriors);
            double bestLogPosterior = (Double.isInfinite(Math.log10(sposteriors[2]))) ? -10000.0 : Math.log10(sposteriors[2]);
            double nextBestLogPosterior = (Double.isInfinite(Math.log10(sposteriors[1]))) ? -10000.0 : Math.log10(sposteriors[1]);
            double lodVsNextBest = bestLogPosterior - nextBestLogPosterior;

            AlleleFrequencyEstimate freq = new AlleleFrequencyEstimate(context.getLocation(),
                                                                       ref,
                                                                       BaseUtils.baseIndexToSimpleBase(altBaseIndex),
                                                                       2,
                                                                       qhat,
                                                                       qstar,
                                                                       (MathUtils.compareDoubles(qhat, reportingBalances[0]) == 0 ? -lodVsNextBest : lodVsRef),
                                                                       lodVsNextBest,
                                                                       pBest,
                                                                       pRef,
                                                                       probs.length,
                                                                       ReadBackedPileup.basePileupAsString(context.getReads(), context.getOffsets()),
                                                                       probs,
                                                                       posteriors);

            return freq;
        }

        return null;
    }

    private double[] getObservationWeights(double[][] probs, int refBaseIndex, int altBaseIndex) {
        return getWeightTableTraces(getWeightTable(probs, refBaseIndex, altBaseIndex, probs.length));
    }

    private double[][] getWeightTable(double[][] probs, int refBaseIndex, int altBaseIndex, int numReadsToConsider) {
        if (numReadsToConsider == 1) {
            double[][] partialProbTable = new double[1][2];
            partialProbTable[0][0] = probs[0][refBaseIndex];
            partialProbTable[0][1] = probs[0][altBaseIndex];

            return partialProbTable;
        }

        double[][] oldPartialProbTable = getWeightTable(probs, refBaseIndex, altBaseIndex, numReadsToConsider - 1);
        double[] traces = getWeightTableTraces(oldPartialProbTable);

        double[][] newPartialProbTable = new double[numReadsToConsider][2];
        for (int row = 0, traceElement = traces.length - 1; row < newPartialProbTable.length; row++, traceElement--) {
            newPartialProbTable[row][0] = traces[traceElement]*probs[numReadsToConsider - 1][refBaseIndex];
            newPartialProbTable[row][1] = traces[traceElement]*probs[numReadsToConsider - 1][altBaseIndex];
        }

        return newPartialProbTable;
    }

    private double[] getWeightTableTraces(double[][] partialProbTable) {
        double[] traces = new double[partialProbTable.length + 1];

        traces[0] = partialProbTable[partialProbTable.length - 1][0];
        traces[partialProbTable.length] = partialProbTable[0][1];

        for (int element = 1; element < traces.length - 1; element++) {
            traces[element] = partialProbTable[partialProbTable.length - element - 1][0] +
                              partialProbTable[partialProbTable.length - element][1];
        }

        return traces;
    }

    private double[] getFractionalCounts(double[][] probs) {
        double[] fractionalCounts = new double[4];

        for (int i = 0; i < probs.length; i++) {
            for (int j = 0; j < 4; j++) {
                fractionalCounts[j] += probs[i][j];
            }
        }

        return fractionalCounts;
    }

    private int getMostFrequentNonRefBase(double[] fractionalCounts, int refBaseIndex) {
        double maxFractionalCounts = -1.0;
        int bestAltBaseIndex = -1;
        for (int altBaseIndex = 0; altBaseIndex < 4; altBaseIndex++) {
            if (altBaseIndex != refBaseIndex) {
                if (fractionalCounts[altBaseIndex] > maxFractionalCounts) {
                    maxFractionalCounts = fractionalCounts[altBaseIndex];
                    bestAltBaseIndex = altBaseIndex;
                }
            }
        }

        return bestAltBaseIndex;
    }

    private AlleleFrequencyEstimate getOneProbAlleleFrequency(char ref, LocusContext context, String rodString) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        ref = Character.toUpperCase(ref);
        
        GenotypeLikelihoods G = new GenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )  {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }
        G.ApplyPrior(ref, Double.NaN);

        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods);
    }

    private String getRodString(RefMetaDataTracker tracker) {
        String rodString = "";

        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )  {
            if ( datum != null )  {
                if ( datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toString();
                } else  {
                    rodString += datum.toSimpleString();
                }
            }
        }
        
        if ( rodString != "" ) { rodString = "[ROD: " + rodString + "]"; }
        return rodString;
    }

    public void setAlleleFrequencyPrior(double freq)
    {
        assert(false);
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(AlleleFrequencyEstimate value, Integer sum) { return 0; }
}
