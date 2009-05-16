
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
    @Argument(fullName="individual_output_prefix", shortName="individual_output_prefix", required=true, doc="prefix to write individual SNP calls to") public String INDIVIDUAL_OUTPUT_PREFIX;

    
    private Random random;
    private SAMFileHeader header;
	private PrintStream discovery_output_file;

    public void initialize() 
    { 
		try
		{
			discovery_output_file = new PrintStream(DISCOVERY_OUTPUT);
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

        random = new Random(42);

		HashSet<String> unique_sample_names = new HashSet<String>();

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
			if (unique_sample_names.contains(sample_name)) { continue; }
			unique_sample_names.add(sample_name);
            sample_names.add(sample_name);
            System.out.println("SAMPLE: " + sample_name);

            SingleSampleGenotyper caller = new SingleSampleGenotyper();
			caller.callsFileName = INDIVIDUAL_OUTPUT_PREFIX + "." + sample_name + ".calls";
			caller.metricsFileName = INDIVIDUAL_OUTPUT_PREFIX + "." + sample_name + ".metrics";
            caller.lodThreshold = lodThreshold;
            caller.fourBaseMode = false;
            caller.printMetrics = false;
            caller.initialize();
            callers.add(caller);
        } 
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {


        // 1. seperate each context.
        LocusContext[] contexts = new LocusContext[sample_names.size()];
        for (int i = 0; i < sample_names.size(); i++)
        {
            contexts[i] = filterLocusContext(context, sample_names.get(i), 0);
        }

        // EM Loop:
	    AlleleFrequencyEstimate[] calls = null;
	    double EM_alt_freq;
	    double EM_N = 0;

		// this line is kinda hacky
		if (MAX_ITERATIONS == 1) { EM_alt_freq = -1; }
		else { EM_alt_freq = 0.5; }

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
	                callers.get(i).setAlleleFrequencyPrior(EM_alt_freq);
	                calls[i] = callers.get(i).map(tracker, ref, contexts[i]);
	                String genotype = calls[i].genotype();

                    likelihood += calls[i].posterior();
	
	
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
	        EM_alt_freq = EM_sum / EM_N;
            trajectory[iterations+1] = EM_alt_freq;
            likelihood_trajectory[iterations+1] = likelihood/(double)EM_N;

            //System.out.printf("DBGTRAJ %f %f %f %f\n", EM_sum, EM_N, trajectory[iterations], trajectory[iterations+1]);
        }

        // 7. Output some statistics.
		double discovery_posterior = 0;
		double discovery_null = 0;
		for (int i = 0; i < sample_names.size(); i++)
		{
			discovery_posterior += calls[i].pBest;
			discovery_null      += calls[i].pRef;
			//System.out.printf("DBG %f %f %c %s\n", calls[i].pBest, calls[i].pRef, ref, calls[i].bases);
		}
		double discovery_lod = discovery_posterior - discovery_null;
		discovery_output_file.printf("%s %f %f %f %f\n", context.getLocation(), EM_alt_freq, discovery_posterior, discovery_null, discovery_lod);

        for (int i = 0; i < sample_names.size(); i++)
        {
            ReadBackedPileup pileup = new ReadBackedPileup(ref, contexts[i]);  
            if (calls[i].depth == 0) { continue; }
			//if (calls[i].lodVsRef < lodThreshold) { continue; }
            out.printf("%s %s %c %f %s %f %f %f %f %f %s\n", context.getLocation(), sample_names.get(i), ref, EM_alt_freq, calls[i].genotype(), calls[i].lodVsRef, calls[i].lodVsNextBest, calls[i].pBest, calls[i].pRef, discovery_lod, pileup.getBases());
        }

        System.out.printf("EVAL %s\n", context.getLocation());
                                
        //for (int i = 0; i < likelihood_trajectory.length; i++)
        //{
        //    System.out.printf("TRAJECTORY %f %f\n", trajectory[i], likelihood_trajectory[i]);
        //}
        //System.out.print("\n\n");

        return null;
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

    private LocusContext filterLocusContext(LocusContext context, String sample_name, int downsample)
    {
        ArrayList<SAMRecord> reads   = new ArrayList<SAMRecord>();
        ArrayList<Integer>   offsets = new ArrayList<Integer>();

        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(header != null);
            //System.out.printf("RG: %s\n", RG);
            assert(header.getReadGroup(RG) != null);

            String sample = header.getReadGroup(RG).getSample();
            if (sample == sample_name) 
            { 
                reads.add(read); 
                offsets.add(offset); 
            }
        }

        if (downsample != 0)
        {
            List<Integer> perm = new ArrayList<Integer>(); 
            for (int i = 0; i < reads.size(); i++) { perm.add(i); }
            perm = Utils.RandomSubset(perm, downsample);
           
            ArrayList<SAMRecord> downsampled_reads = new ArrayList<SAMRecord>();
            ArrayList<Integer> downsampled_offsets = new ArrayList<Integer>();

            for (int i = 0; i < perm.size(); i++)
            {
                downsampled_reads.add(reads.get(perm.get(i)));
                downsampled_offsets.add(offsets.get(perm.get(i)));
            }

            reads = downsampled_reads;
            offsets = downsampled_offsets;
        }

        return new LocusContext(context.getLocation(), reads, offsets);
    }

    public void onTraversalDone() 
    {
        return;
    }

    public String reduceInit() 
    { 
        return "";
    }

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) 
    {
        return "";
    }

    
}
