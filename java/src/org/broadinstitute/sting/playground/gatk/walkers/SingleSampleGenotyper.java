package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods.IndelCall;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import net.sf.samtools.*;

import java.io.*;
import java.util.*;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, String>  
{
    @Argument(fullName="calls",        shortName="calls",        doc="File to write SNP calls to",      required=true) public String callsFileName;
    @Argument(fullName="metrics",      shortName="met",          doc="metrics",      required=false) public String metricsFileName = "/dev/null";
    @Argument(fullName="metInterval",  shortName="mi",           doc="metInterval",  required=false)     public Integer metricsInterval = 50000;
    @Argument(fullName="printMetrics", shortName="printMetrics", doc="printMetrics", required=false)      public Boolean printMetrics = true;
    @Argument(fullName="lodThreshold", shortName="lod",          doc="lodThreshold", required=false)       public Double lodThreshold = 5.0;
    @Argument(fullName="fourBaseMode", shortName="fb",           doc="fourBaseMode", required=false)     public Boolean fourBaseMode = false;
    @Argument(fullName="retest",       shortName="re",           doc="retest",       required=false)     public Boolean retest = false;
    @Argument(fullName="call_indels",  shortName="call_indels",  doc="Call Indels",  required=false)     public Boolean call_indels = false;
    @Argument(fullName="qHom",         shortName="qHom",         doc="qHom",         required=false)      public Double qHom = 0.04;
    @Argument(fullName="qHet",         shortName="qHet",         doc="qHet",         required=false)      public Double qHet = 0.49;
    @Argument(fullName="qHomNonRef",   shortName="qHomNonRef",   doc="qHomNonRef",   required=false)      public Double qHomNonRef = 0.97;
    @Argument(fullName="sample_name_regex", shortName="sample_name_regex", required=false, doc="sample_name_regex") public String SAMPLE_NAME_REGEX = null;

    public AlleleMetrics metrics;
	public PrintStream   calls_file;

	public String sample_name;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) { return true; }
    public boolean requiresReads() { return true; }
    public void initialize() 
	{ 
		try 
		{
			sample_name = null; 
			if (metricsFileName != null) { metrics    = new AlleleMetrics(metricsFileName, lodThreshold); }
			if (callsFileName != null)   { calls_file = new PrintStream(callsFileName); }
		}
		catch (Exception e)
		{
			e.printStackTrace(); 
			System.exit(-1);
		}
	}

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        String rodString = getRodString(tracker);

		if (ref == 'N') { return null; }

		if (context.getReads().size() != 0)
		{
			SAMRecord read = context.getReads().get(0);
	        String RG = (String)(read.getAttribute("RG"));
			SAMReadGroupRecord read_group_record = read.getHeader().getReadGroup(RG);
			if (read_group_record != null)
			{
				String local_sample_name = read.getHeader().getReadGroup(RG).getSample();
				if (SAMPLE_NAME_REGEX != null) { local_sample_name = local_sample_name.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
				if (sample_name == null) { sample_name = local_sample_name; }
				else 
				{ 
					if (! sample_name.equals(local_sample_name)) { System.out.printf("SAMPLE NAME MIXUP: %s vs. %s\n", sample_name, local_sample_name); }
					assert(sample_name.equals(local_sample_name)); 
				}
			}
		}

        AlleleFrequencyEstimate freq = getOneProbAlleleFrequency(ref, context, rodString, sample_name);

        if (printMetrics) {
            if (freq != null) { metrics.nextPosition(freq, tracker); }
            metrics.printMetricsAtLocusIntervals(metricsInterval);
        }

        return freq;
    }

    private AlleleFrequencyEstimate getFourProbAlleleFrequency(char ref, LocusContext context, String rodString, RefMetaDataTracker tracker) {
        double[][] probs = ReadBackedPileup.probDistPileup(context.getReads(), context.getOffsets());
        int refBaseIndex = BaseUtils.simpleBaseToBaseIndex(ref);
        int altBaseIndex = getMostFrequentNonRefBase(getFractionalCounts(probs), refBaseIndex);

        if (refBaseIndex >= 0 && altBaseIndex >= 0) {
            // Set the priors for the hom-ref, het, and non-hom ref (we don't evaluate the
            // possibility of a non-ref het).
            double[] genotypePriors = { 0.999, 1e-3, 1e-5 };

            // I'm not sure what the difference between qhat and qstar is in AlleleMetrics, but
            // to make my life simpler, I evaluate the non-ref balances in genotypeBalances and
            // report the balances in reportingBalances.
            double[] genotypeBalances  = { qHom, qHet, qHomNonRef };
            double[] reportingBalances = { 0.0, 0.5, 1.0 };

            // Compute the probability of the distribution of second-best bases (should be random)
            int[] secondaryBaseCounts = getSecondaryBaseCounts(probs);
            double[] secondaryBaseProbs = new double[4];
            for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                secondaryBaseProbs[baseIndex] = 0.25;
            }
            double secondaryBaseDistProb = MathUtils.multinomialProbability(secondaryBaseCounts, secondaryBaseProbs);

            //System.out.printf("%d %d %d %d ", secondaryBaseCounts[0], secondaryBaseCounts[1], secondaryBaseCounts[2], secondaryBaseCounts[3]);
            //System.out.println(secondaryBaseDistProb);

            // Determine the probability distribution of non-ref alleles
            double[] obsWeights = getObservationWeights(probs, refBaseIndex, altBaseIndex);

            // Convolve them with the binomial sampling probability distribution in order to
            // marginalize the non-ref allele balance.
            double[] posteriors = new double[3];
            double qhat = 0.0, qstar = 0.0, lodVsRef = 0.0, pBest = Double.MIN_VALUE;
            for (int hypothesis = 0; hypothesis < 3; hypothesis++) {
                posteriors[hypothesis] = 0.0;

                for (int weightIndex = 0; weightIndex < obsWeights.length; weightIndex++) {
                    double binomialWeight = MathUtils.binomialProbability(weightIndex, probs.length, genotypeBalances[hypothesis]);

                    // The hom-ref posterior probability distribution may be too weak to
                    // effectively account for noise.  Just set it to equal weighting for
                    // the lower quantile of the non-ref balance range.
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
                                                                       (MathUtils.compareDoubles(qhat, reportingBalances[0]) == 0 ? lodVsNextBest : lodVsRef),
                                                                       lodVsNextBest,
                                                                       pBest,
                                                                       pRef,
                                                                       probs.length,
                                                                       ReadBackedPileup.basePileupAsString(context.getReads(), context.getOffsets()),
                                                                       probs,
                                                                       posteriors,
																	   "unknown_sample");

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
    private int[] getSecondaryBaseCounts(double[][] probs) {
        int[] secondaryBaseCounts = new int[4];

        for (int readIndex = 0; readIndex < probs.length; readIndex++) {
            double bestProb = 0.0;
            for (int baseIndex = 0; baseIndex < probs[readIndex].length; baseIndex++) {
                if (probs[readIndex][baseIndex] > bestProb) {
                    bestProb = probs[readIndex][baseIndex];
                }
            }

            double secondBestProb = 0.0;
            int secondBestBaseIndex = 0;
            for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                if (probs[readIndex][baseIndex] > secondBestProb && probs[readIndex][baseIndex] < bestProb) {
                    secondBestProb = probs[readIndex][baseIndex];
                    secondBestBaseIndex = baseIndex;
                }
            }

            secondaryBaseCounts[secondBestBaseIndex]++;
        }

        return secondaryBaseCounts;
    }

    private AlleleFrequencyEstimate getOneProbAlleleFrequency(char ref, LocusContext context, String rodString, String sample_name) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        ref = Character.toUpperCase(ref);

		// Handle indels.
		if (call_indels)
		{
			String[] indels = BasicPileup.indelPileup(reads, offsets);
			IndelCall indel_call = GenotypeLikelihoods.callIndel(indels);
			if (indel_call != null)
			{
				if (! indel_call.type.equals("ref"))
				{ 
					System.out.printf("INDEL %s %s\n", context.getLocation(), indel_call); 
				}
			}
		}
        
		// Handle single-base polymorphisms.
        GenotypeLikelihoods G = new GenotypeLikelihoods();
        for ( int i = 0; i < reads.size(); i++ )  
        {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }
        G.ApplyPrior(ref, this.allele_frequency_prior);

        if (fourBaseMode) {
            G.applyFourBaseDistributionPrior(pileup.getBases(), pileup.getSecondaryBasePileup());
        }

        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, sample_name);
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

    double allele_frequency_prior = -1;
    public void setAlleleFrequencyPrior(double freq)
    {
        this.allele_frequency_prior = freq;
    }

    // Given result of map function
    private String  confident_ref_interval_contig  = "";
    private long    confident_ref_interval_start   = 0;
    private double  confident_ref_interval_LOD_sum = 0;
    private double  confident_ref_interval_length  = 0;
    private long    last_position_considered       = -1;
    private boolean inside_confident_ref_interval  = false;

    public String reduceInit() 
	{ 
        confident_ref_interval_contig  = "";
        confident_ref_interval_start   = 0;
        confident_ref_interval_LOD_sum = 0;
        confident_ref_interval_length  = 0;
        last_position_considered       = -1;
        inside_confident_ref_interval  = false;
        return "";
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
	{
		calls_file.println(alleleFreq.asTabularString());
		return "";
	}

	/*
    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
	{
        // Print RESULT data for confident calls

       long current_offset = alleleFreq.location.getStart(); //Integer.parseInt(tokens[1]);

        if (inside_confident_ref_interval &&
                ((alleleFreq.lodVsRef > -5.0) || (current_offset != last_position_considered + 1)) )
        {
            // No longer hom-ref, so output a ref line.
            double lod = confident_ref_interval_LOD_sum / confident_ref_interval_length;

            calls_file.format("%s\tCALLER\tREFERENCE\t%d\t%d\t%f\t.\t.\tLENGTH %d\n",
                            confident_ref_interval_contig,
                            confident_ref_interval_start,
                            last_position_considered,
                            lod,
                            (int)(confident_ref_interval_length));

            inside_confident_ref_interval = false;
        }
        else if (inside_confident_ref_interval && (alleleFreq.lodVsRef <= -5.0))
        {
            // Still hom-ref so increment the counters.
            confident_ref_interval_LOD_sum += alleleFreq.lodVsRef;
            confident_ref_interval_length  += 1;
        }
        else if ((!inside_confident_ref_interval) && (alleleFreq.lodVsRef > -5.0))
        {
            // do nothing.
        }
        else if ((!inside_confident_ref_interval) && (alleleFreq.lodVsRef <= -5.0))
        {
            // We moved into a hom-ref region so start a new interval.
            confident_ref_interval_contig  = alleleFreq.location.getContig();
            confident_ref_interval_start   = alleleFreq.location.getStart();
            confident_ref_interval_LOD_sum = alleleFreq.lodVsRef;
            confident_ref_interval_length  = 1;
            inside_confident_ref_interval  = true;
        }

        last_position_considered = current_offset;
        
        if (alleleFreq.lodVsRef >= 5) {
            calls_file.print(alleleFreq.asGFFString());

            //String gtype = genotypeTypeString(alleleFreq.qstar, alleleFreq.N);
            //System.out.print("DEBUG " + gtype + " ");
            //if (gtype.contentEquals("het")) {
            //    System.out.println(alleleFreq.ref + "" + alleleFreq.alt);
            //} else if (gtype.contentEquals("hom")) {
            //    System.out.println(alleleFreq.ref + "" + alleleFreq.ref);
            //} else {
            //    System.out.println(alleleFreq.alt + "" + alleleFreq.alt);
            //}
        }

        return "";
    }
	*/
}
