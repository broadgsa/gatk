
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.util.zip.*;
import java.io.*;

// Draft iterative pooled caller
// j.maguire 4-27-2009

public class PoolCaller extends LocusWalker<AlleleFrequencyEstimate, String> 
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

	AlleleFrequencyEstimate[] calls;
	ArrayList<String>         caller_sums;

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

        /*
        GenomeAnalysisEngine toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getSamReader().getFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();
        */

        sample_names    = new ArrayList<String>();
        callers         = new ArrayList<SingleSampleGenotyper>();
		caller_sums     = new ArrayList<String>();

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
			caller.callsFileName = null;
			caller.metricsFileName = null;
            caller.lodThreshold = lodThreshold;
            caller.fourBaseMode = false;
            caller.printMetrics = false;
			caller.SAMPLE_NAME_REGEX = SAMPLE_NAME_REGEX;
            caller.initialize();
			caller.calls_file = individual_output_file;
			caller_sums.add(caller.reduceInit());
            callers.add(caller);
        } 
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
		if (ref == 'N') { return null; }
		ref = Character.toUpperCase(ref);

        // seperate each context.
        LocusContext forward  = filterLocusContextByStrand(context, "+");
        LocusContext backward = filterLocusContextByStrand(context, "-");

		if (forward.getReads().size() == 0) { return null; }
		if (backward.getReads().size() == 0) { return null; }

		AlleleFrequencyEstimate estimate_both     = EM(tracker, ref, context);
		AlleleFrequencyEstimate estimate_forward  = EM(tracker, ref, forward);
		AlleleFrequencyEstimate estimate_backward = EM(tracker, ref, backward);

		discovery_output_file.printf("%s %c %f %c %f\n",
				estimate_both.asPoolTabularString(),
				estimate_forward.alt,
				estimate_forward.lodVsRef,
				estimate_backward.alt,
				estimate_backward.lodVsRef);
		//discovery_output_file.printf("%s\n", estimate_forward.asPoolTabularString());
		//discovery_output_file.printf("%s\n", estimate_backward.asPoolTabularString());
		//discovery_output_file.printf("\n");

		return null;
	}

	private AlleleFrequencyEstimate EM(RefMetaDataTracker tracker, char ref, LocusContext context)
	{
		if (context.getReads().size() == 0) { return null; }
        LocusContext[] contexts = filterLocusContext(context, sample_names, 0);

        // EM Loop:
	    double EM_alt_freq;
	    double EM_N = 0;
		calls = null;

		// this line is kinda hacky
		if (MAX_ITERATIONS == 1) { EM_alt_freq = -1; }
		else { EM_alt_freq = 0.5; }

        // (this loop is the EM cycle)
        double[] trajectory = new double[MAX_ITERATIONS + 1]; trajectory[0] = EM_alt_freq;
        double[] likelihood_trajectory = new double[MAX_ITERATIONS + 1]; likelihood_trajectory[0] = 0.0;
        boolean is_a_snp = false;

		// Pick the alternate base
		char alt = 'N';
		{
        	ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
			String bases = pileup.getBases();
			int A = 0;
			int C = 0;
			int G = 0; 
			int T = 0;
			int max_count = -1;
			for (int i = 0; i < bases.length(); i++)
			{
				char b = bases.charAt(i);
				if (b == ref) { continue; }
				switch (b)
				{
					case 'A' : A += 1; if (A > max_count) { max_count = A; alt = 'A'; } break;
					case 'C' : C += 1; if (C > max_count) { max_count = C; alt = 'C'; } break;
					case 'G' : G += 1; if (G > max_count) { max_count = G; alt = 'G'; } break;
					case 'T' : T += 1; if (T > max_count) { max_count = T; alt = 'T'; } break;
				}
			}
		}

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
	                		//System.out.printf("DBG: %s %f %f\n",
							//				context.getLocation(),
	                        //              calls[i].lodVsNextBest,
	                        //              calls[i].lodVsRef);
				            EM_sum += calls[i].emperical_allele_frequency() * calls[i].N;
				            EM_N   += calls[i].N;
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

            //System.out.printf("DBGTRAJ %s %f %f %f %f %f %f\n", 
			//						context.getLocation(),
			//						EM_sum, 
			//						EM_N, 
			//						trajectory[iterations], 
			//						trajectory[iterations+1], 
			//						likelihood_trajectory[iterations],
			//						likelihood_trajectory[iterations+1]);
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

			discovery_likelihood += calls[i].pBest;
			discovery_null       += calls[i].pRef;
			//System.out.printf("DBG %f %f %c %s\n", calls[i].pBest, calls[i].pRef, ref, calls[i].bases);
		
			if (calls[i].qhat == 0.0) { n_ref += 1; }	
			if (calls[i].qhat == 0.5) { n_het += 1; }	
			if (calls[i].qhat == 1.0) { n_hom += 1; }	
		}
		double discovery_lod = (discovery_likelihood + discovery_prior) - (discovery_null + discovery_null_prior);
		if (discovery_lod <= 0) { alt = 'N'; }
		//discovery_output_file.printf("%s %c %c %f %f %f %f %f %f %d %d %d\n", context.getLocation(), ref, alt, EM_alt_freq, discovery_likelihood, discovery_null, discovery_prior, discovery_lod, EM_N, n_ref, n_het, n_hom);

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
		return estimate;


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

    private void CollectStrandInformation(char ref, LocusContext context)
	{
		List<SAMRecord> reads   = context.getReads();
		List<Integer>   offsets = context.getOffsets();

	}

    public void onTraversalDone(String result) 
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
		return;
    }

    public String reduceInit() 
    { 
		discovery_output_file.printf("loc ref alt EM_alt_freq discovery_likelihood discovery_null discovery_prior discovery_lod EM_N n_ref n_het n_hom fw_alt fw_lod bw_alt bw_lod\n");
		for (int i = 0; i < callers.size(); i++)
		{
			callers.get(i).reduceInit(); 
		}
		return "";
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
    {
		if (calls == null) { return ""; }
		for (int i = 0; i < callers.size(); i++)
		{
			if (calls == null) { System.err.printf("calls == null\n"); }
			if (calls[i] == null) { System.err.printf("calls[%d] == null\n", i); }
			if (caller_sums == null) { System.err.printf("caller_sums == null\n"); }
			if (callers.get(i) == null) { System.err.printf("callers[%d] == null\n", i); }
			if (caller_sums.get(i) == null) { System.err.printf("caller_sums[%d] == null\n", i); }
			caller_sums.set(i, callers.get(i).reduce(calls[i], caller_sums.get(i))); 
		}
        return "";
    }

    
}
