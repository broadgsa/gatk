package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;

import net.sf.samtools.SAMRecord;

import java.util.*;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, Integer>  {
    @Argument(fullName="fourBaseMode",required=false,defaultValue="false")
    public Boolean fourBaseMode;

    @Argument(fullName="decideOnBase",required=false,defaultValue="false")
    public Boolean decideOnBase;
    
    private AlleleMetrics metrics;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) { return true; } // We are keeping all the reads
    public boolean requiresReads() { return true; }
    public void initialize() { metrics = new AlleleMetrics("metrics.out"); }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        String rodString = getRodString(tracker);

        AlleleFrequencyEstimate freq = null;
        if (fourBaseMode) 
        {
            // Compute four-base prob genotype likelihoods
            freq = getFourProbAlleleFrequency(ref, context, rodString);
        } 
        else if (decideOnBase) 
        {
        } 
        else 
        {
            // Compute single quality score genotype likelihoods
            freq = getOneProbAlleleFrequency(ref, context, rodString);
        }

        if (freq != null) { metrics.nextPosition(freq, tracker); }
        metrics.printMetricsAtLocusIntervals(1000);

        return freq;
    }

    private AlleleFrequencyEstimate getFourProbAlleleFrequency(char ref, LocusContext context, String rodString) {
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
            System.out.println(context.getLocation().toString() + " " + refBaseIndex + " " + altBaseIndex + " " + probs.length + " " + rodString);
            for (int i = 0; i < probs.length; i++) {
                System.out.printf("  [ %4.4f %4.4f %4.4f %4.4f ]\n", probs[i][0], probs[i][1], probs[i][2], probs[i][3]);
            }
            
            double[] obsWeights = getObservationWeights(probs, refBaseIndex, altBaseIndex);

            System.out.print("  Weights: ");
            for (int i = 0; i < obsWeights.length; i++) {
                System.out.printf("%4.4f ", obsWeights[i]);
            }
            System.out.println();

            double[] genotypePriors   = { 0.999, 1e-3, 1e-5 };
            double[] genotypeBalances = { 0.02,  0.5,  0.98 };

            double[] posteriors = new double[3];
            double qhat = 0.0, qstar = 0.0, lodVsRef = 0.0, lodVsNextBest = 0.0, pBest = Double.MIN_VALUE;

            for (int hypothesis = 0; hypothesis < 3; hypothesis++) {
                posteriors[hypothesis] = 0.0;

                for (int weightIndex = 0; weightIndex < obsWeights.length; weightIndex++) {
                    posteriors[hypothesis] += obsWeights[weightIndex]*binomialProb(weightIndex, probs.length, genotypeBalances[hypothesis]);
                }
                posteriors[hypothesis] *= genotypePriors[hypothesis];

                System.out.printf("  Hypothesis %d %f %f %f\n", hypothesis, genotypeBalances[hypothesis], posteriors[hypothesis], Math.log10(posteriors[hypothesis]/posteriors[0]));

                if (posteriors[hypothesis] > pBest) {
                    qhat = genotypeBalances[hypothesis];
                    qstar = genotypeBalances[hypothesis];
                    lodVsRef = Math.log10(posteriors[hypothesis]/posteriors[0]);
                    pBest = posteriors[hypothesis];
                }
            }

            double pRef = posteriors[0];
            
            System.out.println("\n");

            return new AlleleFrequencyEstimate(context.getLocation(),
                                               ref,
                                               BaseUtils.baseIndexToSimpleBase(altBaseIndex),
                                               2,
                                               qhat,
                                               qstar,
                                               lodVsRef,
                                               lodVsNextBest,
                                               pBest,
                                               pRef,
                                               probs.length,
                                               ReadBackedPileup.basePileupAsString(context.getReads(), context.getOffsets()),
                                               probs,
                                               posteriors);

        }

        return null;
    }

    double binomialProb(long k, long n, double p) {
        // k - number of successes
        // n - number of Bernoulli trials
        // p - probability of success

        if ((n*p < 5) && (n*(1-p) < 5))
        {
            // For small n and the edges, compute it directly.
            return (double)nchoosek(n, k) * Math.pow(p, k) * Math.pow(1-p, n-k);
        }
        else
        {
            // For large n, approximate with a gaussian.
            double mean  = (double)(n*p);
            double var   = Math.sqrt((double)(n*p)*(1.0-p));
            double ans   = (double)(1.0 / (var*Math.sqrt(2*Math.PI)))*Math.exp(-1.0 * Math.pow((double)k-mean,2)/(2.0*var*var));
            double check = (double)nchoosek(n, k) * Math.pow(p, k) * Math.pow(1-p, n-k);
            double residual = ans - check;

            return check;
        }
    }

    long nchoosek(long n, long k) {
        long m = n - k;
        if (k < m)
            k = m;

        long t = 1;
        for (long i = n, j = 1; i > k; i--, j++)
            t = t * i / j;
        return t;
    }

    private double[] getObservationWeights(double[][] probs, int refBaseIndex, int altBaseIndex) {
        //if (probs.length <= 10) 
        //{
        //    return getWeightTableTraces(getWeightTable(probs, refBaseIndex, altBaseIndex, probs.length));
        //}
        //else
        //{
            return FastObservationWeights(probs, refBaseIndex, altBaseIndex);
        //}
    }

    private double[] FastObservationWeights(double[][] probs, int ref, int alt)
    {
        List<int[]> paths        = new ArrayList<int[]>();
        List<Double> likelihoods = new ArrayList<Double>();

        List<int[]>  output_paths       = new ArrayList<int[]>();
        List<Double> output_likelihoods = new ArrayList<Double>();

        HashMap<String,Boolean> done_paths = new HashMap<String,Boolean>();

        // 1. Find the best path.
        int[] best_path = new int[probs.length];
        double best_likelihood = 0;
        for (int i = 0; i < probs.length; i++)
        {
            int max;
            double max_p;
            if (probs[i][ref] >= probs[i][alt]) { max = ref; max_p = probs[i][ref]; }
            else                                { max = alt; max_p = probs[i][alt]; }
            best_path[i] = max;
            best_likelihood += Math.log10(max_p);
        }
        output_paths.add(best_path);
        output_likelihoods.add(best_likelihood);
        String s = ""; for (int j = 0; j < best_path.length; j++) { s += best_path[j]; }
        done_paths.put(s,true);

        // 2. Enumerate all paths one-away from the best path
        for (int i = 0; i < best_path.length; i++)
        {
            int[] path = Arrays.copyOf(best_path, best_path.length);
            double likelihood;
            if (path[i] == ref)
            {
                path[i] = alt;
                likelihood = best_likelihood - Math.log10(probs[i][ref]) + Math.log10(probs[i][alt]);
            }
            else
            {
                path[i] = ref;
                likelihood = best_likelihood - Math.log10(probs[i][alt]) + Math.log10(probs[i][ref]);
            }
            paths.add(path);
            likelihoods.add(likelihood);
        }

        // 3. Sort paths by likelihood
        Integer[] permutation = Utils.SortPermutation(likelihoods);
        paths = Utils.PermuteList(paths, permutation);
        likelihoods = Utils.PermuteList(likelihoods, permutation);

        while ((output_paths.size() < 10) && (paths.size() > 0))
        {
	        // 4. Choose the next best path
            int[] next_best_path = paths.get(paths.size()-1);
            double next_best_likelihood = likelihoods.get(likelihoods.size()-1);
	        output_paths.add(next_best_path);
	        output_likelihoods.add(next_best_likelihood);
            paths.remove(paths.size()-1);
            likelihoods.remove(likelihoods.size()-1);
            s = ""; for (int j = 0; j < next_best_path.length; j++) { s += next_best_path[j]; }
            done_paths.put(s,true);

                /*
            if (likelihoods.get(likelihoods.size()-1) < probs.length*Math.log10(0.5)) 
            {
                break;
            }
                */

	        // 5. Enumerate all paths one-away from next best
	        for (int i = 0; i < best_path.length; i++)
	        {
	            int[] path = Arrays.copyOf(next_best_path, next_best_path.length);
	            double likelihood;
	            if (path[i] == ref)
	            {
	                path[i] = alt;
	                likelihood = next_best_likelihood - Math.log10(probs[i][ref]) + Math.log10(probs[i][alt]);
	            }
	            else
	            {
	                path[i] = ref;
	                likelihood = next_best_likelihood - Math.log10(probs[i][alt]) + Math.log10(probs[i][ref]);
	            }

                s = "";
                for (int j = 0; j < path.length; j++) { s += path[j]; }
                if (done_paths.get(s) == null) 
                {
	                paths.add(path);
	                likelihoods.add(likelihood);
                }
	        }
	
	        // 6. Re-sort
	        permutation = Utils.SortPermutation(likelihoods);
            paths = Utils.PermuteList(paths, permutation);
            likelihoods = Utils.PermuteList(likelihoods, permutation);
        }

        for (int i = 0; i < output_paths.size(); i++)
        {
            int[] path = output_paths.get(i);
            double likelihood = output_likelihoods.get(i);

            System.out.printf("DBG %d ", i);
            for (int j = 0; j < path.length; j++)
            {
                System.out.printf("%c", BaseUtils.baseIndexToSimpleBase(path[j]));
            }
            System.out.printf(" %f\n", likelihood);
        }
        System.out.printf("\n");

        double[] ans = new double[probs.length+1];
        for (int i = 0; i < ans.length; i++) { ans[i] = 1.0; }
        return ans;
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

        System.out.printf("%s %s %s %s\n", context.getLocation(), ref, bases, G.toString(ref), rodString);

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
