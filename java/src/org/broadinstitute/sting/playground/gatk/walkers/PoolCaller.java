
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.util.zip.*;
import java.io.*;

// Draft iterative pooled caller
// j.maguire 4-27-2009

public class PoolCaller extends LocusWalker<AlleleFrequencyEstimate[], String[]> 
{
    List<SingleSampleGenotyper> callers = null;
    List<String> sample_names = null;

    @Argument(required=false, shortName="fractional_counts", doc="should we use fractional counts?") public boolean FRACTIONAL_COUNTS = false;
    @Argument(required=false, shortName="max_iterations", doc="Maximum number of iterations for EM") public int MAX_ITERATIONS = 10;
    @Argument(fullName="lodThreshold", shortName="lod", required=false, doc="lod threshold for outputting individual genotypes")       public Double lodThreshold = 2.0;
    @Argument(fullName="discovery_output", shortName="discovery_output", required=true, doc="file to write SNP discovery output to")       public String DISCOVERY_OUTPUT;
    @Argument(fullName="individual_output", shortName="individual_output", required=true, doc="file to write individual SNP calls to") public String INDIVIDUAL_OUTPUT;
    @Argument(fullName="sample_name_regex", shortName="sample_name_regex", required=false, doc="sample_name_regex") public String SAMPLE_NAME_REGEX = null;

    private Random random;
    private SAMFileHeader header;
	private PrintStream discovery_output_file;
	private PrintStream individual_output_file;

    public void initialize() 
    { 
		try
		{
			discovery_output_file = new PrintStream(DISCOVERY_OUTPUT);
			individual_output_file = new PrintStream(new GZIPOutputStream(new FileOutputStream(INDIVIDUAL_OUTPUT)));
			individual_output_file.println(AlleleFrequencyEstimate.asTabularStringHeader());
		} 
		catch (Exception e)
		{
			e.printStackTrace(); 
			System.exit(-1);
		}

        GenomeAnalysisEngine toolkit = this.getToolkit();
        this.header = toolkit.getEngine().getSAMHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        sample_names    = new ArrayList<String>();
        callers         = new ArrayList<SingleSampleGenotyper>();

        random = new Random(42);

		HashSet<String> unique_sample_names = new HashSet<String>();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();

			if (SAMPLE_NAME_REGEX != null) { sample_name = sample_name.replaceAll(SAMPLE_NAME_REGEX, "$1"); }

			if (unique_sample_names.contains(sample_name)) { continue; }
			unique_sample_names.add(sample_name);
            sample_names.add(sample_name);
            System.out.println("SAMPLE: " + sample_name);

            SingleSampleGenotyper caller = new SingleSampleGenotyper();
			caller.VARIANTS_FILE = null;
			caller.METRICS_FILE = null;
            caller.LOD_THRESHOLD = lodThreshold;
            caller.IGNORE_SECONDARY_BASES = true;
            caller.SUPPRESS_METRICS = true;
			caller.SAMPLE_NAME_REGEX = SAMPLE_NAME_REGEX;
            caller.initialize();
			caller.variantsOut = individual_output_file;
            callers.add(caller);
        } 
    }

    public AlleleFrequencyEstimate[] map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
		if (ref == 'N') { return null; }
		ref = Character.toUpperCase(ref);

        // seperate each context.
        LocusContext forward  = filterLocusContextByStrand(context, "+");
        LocusContext backward = filterLocusContextByStrand(context, "-");

		if (forward.getReads().size() == 0) { return null; }
		if (backward.getReads().size() == 0) { return null; }

		// Pick the alternate base
		char alt = 'N';
		{
			EM_Result result_both = EM(tracker, ref, context, -1, 'N', 1, lodThreshold, callers);
			int[] counts = new int[4];
			if (result_both.individuals == null) { return null; }
			for (int i = 0; i < result_both.individuals.length; i++)
			{
				if (result_both.individuals[i] == null) { continue; }
				if (result_both.individuals[i].lodVsRef >= lodThreshold) 
				{
					counts[BaseUtils.simpleBaseToBaseIndex(result_both.individuals[i].alt)] += 1;
				}
				Integer[] perm = Utils.SortPermutation(counts);
				alt = BaseUtils.baseIndexToSimpleBase(perm[3]);	
			}
		}

		double EM_alt_freq;
		if (MAX_ITERATIONS == 1) { EM_alt_freq = -1; }
		else                     { EM_alt_freq = 0.5; }

		EM_Result result_both     = EM(tracker, ref, context, EM_alt_freq, alt, MAX_ITERATIONS, lodThreshold, callers);
		EM_Result result_forward  = EM(tracker, ref, forward, EM_alt_freq, alt, MAX_ITERATIONS, lodThreshold, callers);
		EM_Result result_backward = EM(tracker, ref, backward, EM_alt_freq, alt, MAX_ITERATIONS, lodThreshold, callers);

		EM_Result null_both       = EM(tracker, ref, context,  0, alt, 1, 1e-3, callers);
		EM_Result null_forward    = EM(tracker, ref, forward,  0, alt, 1, 1e-3, callers);
		EM_Result null_backward   = EM(tracker, ref, backward, 0, alt, 1, 1e-3, callers);

		if (result_both.pool == null) { return null; }
		AlleleFrequencyEstimate estimate_both = result_both.pool;

		double lod_forward;
		double lod_backward;
		double lod_both;
		double strand_score;
		char forward_alt;
		char backward_alt;

		if ((result_forward.pool == null) ||
		    (result_backward.pool == null) ||
			(null_both == null) ||
			(null_forward == null) ||
		    (null_backward == null))
		{
			lod_forward = 0;
			lod_backward = 0;
			lod_both = 0;
			strand_score = 0;
			forward_alt = 'N';
			backward_alt = 'N';
		}
		else
		{
			AlleleFrequencyEstimate estimate_forward       = result_forward.pool;
			AlleleFrequencyEstimate estimate_backward      = result_backward.pool;

	
			double p_D_both     = 0;
			double p_D_forward  = 0;
			double p_D_backward = 0;
			double p_D_null_both = 0;
			double p_D_null_forward = 0;
			double p_D_null_backward = 0;
			for (int i = 0; i < result_both.individuals.length; i++)
			{
				double sum_both = 0;
				double sum_forward = 0;
				double sum_backward = 0;
	
				double sum_null_both     = 0;
				double sum_null_forward  = 0;
				double sum_null_backward = 0;
	
				for (int j = 0; j < result_both.individuals[i].genotypeLikelihoods.likelihoods.length; j++)
				{
					sum_both     += Math.pow(10, result_both.individuals[i].genotypeLikelihoods.likelihoods[j]);
					sum_forward  += Math.pow(10, result_forward.individuals[i].genotypeLikelihoods.likelihoods[j]);
					sum_backward += Math.pow(10, result_backward.individuals[i].genotypeLikelihoods.likelihoods[j]);
					sum_null_both += Math.pow(10, null_both.individuals[i].genotypeLikelihoods.likelihoods[j]);
					sum_null_forward += Math.pow(10, null_forward.individuals[i].genotypeLikelihoods.likelihoods[j]);
					sum_null_backward += Math.pow(10, null_backward.individuals[i].genotypeLikelihoods.likelihoods[j]);
				}
	
				p_D_both     += Math.log10(sum_both);
				p_D_forward  += Math.log10(sum_forward);
				p_D_backward += Math.log10(sum_backward);
	
				p_D_null_both     += Math.log10(sum_null_both);
				p_D_null_forward  += Math.log10(sum_null_forward);
				p_D_null_backward += Math.log10(sum_null_backward);
			}
			forward_alt = estimate_forward.alt;
			backward_alt = estimate_backward.alt;
			lod_forward  = (p_D_forward  + p_D_null_backward) - p_D_null_both;
			lod_backward = (p_D_backward + p_D_null_backward) - p_D_null_both;
			lod_both     = p_D_both - p_D_null_both;
			strand_score = Math.max(lod_forward - lod_both, lod_backward - lod_both);
		}

		System.out.printf("DBG %s %f %f %f %f\n", context.getLocation(), result_both.pool.pBest, null_both.pool.pBest, result_both.pool.pRef, null_both.pool.pBest);

		discovery_output_file.printf("%s %c %c %f\n",
				estimate_both.asPoolTabularString(),
				forward_alt,
				backward_alt,
				strand_score);

		return result_both.individuals;
	}

	private class EM_Result
	{
		AlleleFrequencyEstimate pool;
		AlleleFrequencyEstimate[] individuals;

		public EM_Result(AlleleFrequencyEstimate pool, AlleleFrequencyEstimate[] individuals) 
		{
			this.pool = pool;
			this.individuals = individuals;
		}

		// Construct an EM_Result that indicates no data.
		public EM_Result()
		{
			this.pool = null; 
			this.individuals = null;
		}
	}

	private EM_Result EM(RefMetaDataTracker tracker, char ref, LocusContext context, double EM_alt_freq, char alt, int MAX_ITERATIONS, double lodThreshold, List<SingleSampleGenotyper> callers)
	{
		if (context.getReads().size() == 0) { return null; }
        LocusContext[] contexts = filterLocusContext(context, sample_names, 0);

        // EM Loop:
	    double EM_N = 0;
		AlleleFrequencyEstimate[] calls = null;

        // (this loop is the EM cycle)
        double[] trajectory = new double[MAX_ITERATIONS + 1]; trajectory[0] = EM_alt_freq;
        double[] likelihood_trajectory = new double[MAX_ITERATIONS + 1]; likelihood_trajectory[0] = 0.0;
        boolean is_a_snp = false;

        for (int iterations = 0; iterations < MAX_ITERATIONS; iterations++)
        { 
	        // 6. Re-call from shallow coverage using the estimated frequency as a prior, 
	        //    and compare to true deep calls, 
	        //    and compute new MAF estimate.
	        calls = new AlleleFrequencyEstimate[sample_names.size()];
	        EM_N        = 0.0;
            double EM_sum = 0.0;
            double likelihood = 0.0;
            is_a_snp = false;

	        for (int i = 0; i < sample_names.size(); i++)
	        {
	                callers.get(i).setAlleleFrequencyPrior(EM_alt_freq, alt);
	                calls[i] = callers.get(i).map(tracker, ref, contexts[i]);
	                String genotype = calls[i].genotype();

					if (calls[i].depth == 0) { continue; }

                    likelihood += calls[i].pBest;
	
                    if (! FRACTIONAL_COUNTS)
                    {
							if (Math.abs(calls[i].lodVsRef) >= lodThreshold) 
							{
								EM_sum += calls[i].emperical_allele_frequency() * calls[i].N;
					            EM_N   += calls[i].N;
							}
                    }
                    else
                    {
	                    for (int j = 0; j <= calls[i].N; j++)
	                    {
                            if (Double.isInfinite(calls[i].posteriors[j])) { calls[i].posteriors[j] = -10000; }
	                        System.out.printf("DBG3: %d %f %d\n", j, calls[i].posteriors[j], calls[i].N);
	                        EM_sum += Math.pow(10,calls[i].posteriors[j]) * (double)j;
	                        EM_N   += calls[i].N;
	                    }
                    }   

            }
	        EM_alt_freq = Math.min(EM_sum / EM_N, 0.99999);
            trajectory[iterations+1] = EM_alt_freq;
            likelihood_trajectory[iterations+1] = likelihood;

			if (likelihood_trajectory[iterations] == likelihood_trajectory[iterations+1]) { break; }

			/*
            System.out.printf("DBGTRAJ %s %f %f %f %f %f %f\n", 
									context.getLocation(),
									EM_sum, 
									EM_N, 
									trajectory[iterations], 
									trajectory[iterations+1], 
									likelihood_trajectory[iterations],
									likelihood_trajectory[iterations+1]);
			*/
        }

        // 7. Output some statistics.
		double discovery_likelihood = 0;
		double discovery_null = 0;
		double discovery_prior = 0;
		double discovery_null_prior = 0;
		int n_ref = 0;
		int n_het = 0;
		int n_hom = 0;
		for (int i = 1; i < EM_N-1; i++)
		{
			//discovery_prior += 1e-3/(double)i;
		}
		//discovery_prior = Math.log10(discovery_prior);
		//discovery_null_prior = Math.log10(1.0 - discovery_prior);
		for (int i = 0; i < sample_names.size(); i++)
		{
			if (calls[i].depth == 0) { continue; }

			if (calls[i].lodVsRef < lodThreshold) { continue; }

			discovery_likelihood += calls[i].pBest;
			discovery_null       += calls[i].pRef;
		
			if (calls[i].qhat == 0.0) { n_ref += 1; }	
			if (calls[i].qhat == 0.5) { n_het += 1; }	
			if (calls[i].qhat == 1.0) { n_hom += 1; }	
		}
		double discovery_lod = (discovery_likelihood + discovery_prior) - (discovery_null + discovery_null_prior);
		if (discovery_lod <= 0) { alt = 'N'; }
		//discovery_output_file.printf("%s %c %c %f %f %f %f %f %f %d %d %d\n", context.getLocation(), ref, alt, EM_alt_freq, discovery_likelihood, discovery_null, discovery_prior, discovery_lod, EM_N, n_ref, n_het, n_hom);

		if (EM_N == 0) { return new EM_Result(); }

		AlleleFrequencyEstimate estimate = new AlleleFrequencyEstimate(context.getLocation(), 
																			ref, 
																			alt, 
																			(int)EM_N, 
																			EM_alt_freq, 
																			EM_alt_freq, 
																			discovery_lod,
																			0.0,
																			discovery_likelihood + discovery_prior,
																			discovery_null + discovery_prior,
																			context.getReads().size(),
																			(String)null,
																			(double[][])null,
																			(double[])null,
																			(String)null);
		estimate.n_ref = n_ref; // HACK
		estimate.n_het = n_het; // HACK
		estimate.n_hom = n_hom; // HACK
		return new EM_Result(estimate, calls);


        //for (int i = 0; i < likelihood_trajectory.length; i++)
        //{
        //    System.out.printf("TRAJECTORY %f %f\n", trajectory[i], likelihood_trajectory[i]);
        //}
        //System.out.print("\n\n");
    }

    private LocusContext poolLocusContext(LocusContext[] contexts)
    {
        ArrayList<SAMRecord> reads   = new ArrayList<SAMRecord>();
        ArrayList<Integer>   offsets = new ArrayList<Integer>();

        GenomeLoc location = null;

        for (int i = 0; i < contexts.length; i++)
        {
            if (contexts[i] != null)
            {
                location = contexts[i].getLocation();
	            reads.addAll(contexts[i].getReads());
	            offsets.addAll(contexts[i].getOffsets());
            }
        }

        return new LocusContext(location, reads, offsets);
    }

	private LocusContext filterLocusContextByStrand(LocusContext context, String strand)
	{
		ArrayList<SAMRecord> reads = new ArrayList<SAMRecord>();
		ArrayList<Integer> offsets = new ArrayList<Integer>();

		for (int i = 0; i < context.getReads().size(); i++)
		{
			SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);

			// Filter for strandedness
			if ((!strand.contains("+")) && (read.getReadNegativeStrandFlag() == false)) { continue; }
			if ((!strand.contains("-")) && (read.getReadNegativeStrandFlag() == true))  { continue; }
			reads.add(read);
			offsets.add(offset);
		}
		return new LocusContext(context.getLocation(), reads, offsets);
	}

    private LocusContext[] filterLocusContext(LocusContext context, List<String> sample_names, int downsample)
    {
		HashMap<String,Integer> index = new HashMap<String,Integer>();
		for (int i = 0; i < sample_names.size(); i++)
		{
			index.put(sample_names.get(i), i);
		}

		LocusContext[] contexts = new LocusContext[sample_names.size()];
		ArrayList<SAMRecord>[] reads = new ArrayList[sample_names.size()];
		ArrayList<Integer>[] offsets = new ArrayList[sample_names.size()];

		for (int i = 0; i < sample_names.size(); i++)
		{
			reads[i] = new ArrayList<SAMRecord>();
			offsets[i] = new ArrayList<Integer>();
		}

        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(header != null);
            assert(header.getReadGroup(RG) != null);

            String sample = header.getReadGroup(RG).getSample();
			if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }
            reads[index.get(sample)].add(read); 
            offsets[index.get(sample)].add(offset); 
        }

        if (downsample != 0)
        {
			for (int j = 0; j < reads.length; j++)
			{
	            List<Integer> perm = new ArrayList<Integer>(); 
	            for (int i = 0; i < reads[j].size(); i++) { perm.add(i); }
	            perm = Utils.RandomSubset(perm, downsample);
	           
	            ArrayList<SAMRecord> downsampled_reads = new ArrayList<SAMRecord>();
	            ArrayList<Integer> downsampled_offsets = new ArrayList<Integer>();
	
	            for (int i = 0; i < perm.size(); i++)
	            {
	                downsampled_reads.add(reads[j].get(perm.get(i)));
	                downsampled_offsets.add(offsets[j].get(perm.get(i)));
	            }
	
	            reads[j] = downsampled_reads;
	            offsets[j] = downsampled_offsets;
				contexts[j] = new LocusContext(context.getLocation(), reads[j], offsets[j]);
			}
        }
		else
		{
			for (int j = 0; j < reads.length; j++)
			{
				contexts[j] = new LocusContext(context.getLocation(), reads[j], offsets[j]);
			}
		}

        return contexts;
    }

    public void onTraversalDone(String[] result) 
    {
		try 
		{
			discovery_output_file.flush();
			discovery_output_file.close();
			individual_output_file.flush();
			individual_output_file.close();
		}
		catch (Exception e)
		{
			e.printStackTrace(); 
		}
		out.println("PoolCaller done.\n");
		return;
    }

	public String[] single_sample_reduce_sums = null;
    public String[] reduceInit() 
    { 
		discovery_output_file.printf("loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom fw_alt bw_alt strand_score\n");
		String[] single_sample_reduce_sums = new String[callers.size()];
		for (int i = 0; i < callers.size(); i++)
		{
			single_sample_reduce_sums[i] = callers.get(i).reduceInit(); 
		}
		return single_sample_reduce_sums;
    }

    public String[] reduce(AlleleFrequencyEstimate[] alleleFreqs, String[] sum) 
    {
		if (alleleFreqs == null) { return sum; }
		for (int i = 0; i < callers.size(); i++)
		{
			sum[i] = callers.get(i).reduce(alleleFreqs[i], sum[i]);
		}
        return sum;
    }

    
}
