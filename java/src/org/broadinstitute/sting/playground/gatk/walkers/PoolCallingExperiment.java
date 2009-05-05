
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;

public class PoolCallingExperiment extends LocusWalker<AlleleFrequencyEstimate, String> 
{
    List<AlleleFrequencyWalker> deep_callers = null;
    List<AlleleFrequencyWalker> shallow_callers = null;
    AlleleFrequencyWalker       pooled_caller = null;
    List<String> sample_names = null;

    @Argument(required=false, shortName="downsample", doc="downsample", defaultValue="4") public int DOWNSAMPLE;
    @Argument(required=false, shortName="downsample_noise", doc="downsample noise", defaultValue="3") public int DOWNSAMPLE_NOISE;
    @Argument(required=false, shortName="log_metrics", doc="log metrics", defaultValue="true") public boolean LOG_METRICS;
    @Argument(required=false, shortName="fractional_counts", doc="fractional counts", defaultValue="false") public boolean FRACTIONAL_COUNTS;
    
    private Random random;

    public void initialize() 
    { 
        GenomeAnalysisTK toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getEngine().getSAMHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        sample_names    = new ArrayList<String>();
        deep_callers    = new ArrayList<AlleleFrequencyWalker>();
        shallow_callers = new ArrayList<AlleleFrequencyWalker>();

        random = new Random(42);

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            //System.out.println("SAMPLE: " + sample_name);

            AlleleFrequencyWalker deep    = new AlleleFrequencyWalker();
            AlleleFrequencyWalker shallow = new AlleleFrequencyWalker();

            deep.N = 2;
            deep.DOWNSAMPLE = 0;
            deep.GFF_OUTPUT_FILE = "/dev/null";
            deep.initalize();

            shallow.N = 2;
            shallow.DOWNSAMPLE = 0;
            shallow.GFF_OUTPUT_FILE = "/dev/null";
            shallow.initalize();

            sample_names.add(sample_name);
            deep_callers.add(deep);
            shallow_callers.add(shallow);
        } 

        pooled_caller = new AlleleFrequencyWalker();
        pooled_caller.N = sample_names.size() * 2;
        pooled_caller.DOWNSAMPLE = 0;
        pooled_caller.GFF_OUTPUT_FILE = "/dev/null";
        pooled_caller.initalize();
        pooled_caller.reduceInit();
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        // 1. seperate each context.
        LocusContext[] deep_contexts = new LocusContext[sample_names.size()];
        for (int i = 0; i < sample_names.size(); i++)
        {
            deep_contexts[i] = filterLocusContext(context, sample_names.get(i), 0);
        }

        // 2. Pool the deep contexts. (for the hell of it)
        LocusContext deep_pooled_context = poolLocusContext(deep_contexts);
        
        // 3. Call individuals from deep coverage.
        AlleleFrequencyEstimate[] deep_calls    = new AlleleFrequencyEstimate[sample_names.size()];
        for (int i = 0; i < sample_names.size(); i++)
        {
            deep_calls[i] = deep_callers.get(i).map(tracker, ref, deep_contexts[i]);
        }

        // 4. Count "true" allele frequency from deep calls, and compare to the estimated frequency from the pool.
        // 4.1. Count "true" allele frequency from deep calls, and downsample the individuals with solid truth.
        double true_alt_freq = 0.0;
        double true_N        = 0.0;
        LocusContext[] downsampled_contexts = new LocusContext[sample_names.size()];
        for (int i = 0; i < deep_calls.length; i++)
        {
            downsampled_contexts[i] = null;
            if ((deep_calls[i].lodVsNextBest >= 5.0) || (deep_calls[i].lodVsRef <= -5.0))
            {
                true_alt_freq += deep_calls[i].emperical_allele_frequency() * deep_calls[i].N;
                true_N        += 2;
                downsampled_contexts[i] = filterLocusContext(context, sample_names.get(i), DOWNSAMPLE + (int)(random.nextGaussian()*DOWNSAMPLE_NOISE));
            }
        }
        true_alt_freq /= true_N;

        if (true_N == 0.0) { return null; } // just bail if true_N == 0.

        // 4.2. Pool just the contexts that have truth data.
        LocusContext pooled_context = poolLocusContext(downsampled_contexts);


        // EM Loop:
        AlleleFrequencyEstimate pooled_call = null;

	    double correct_shallow_calls = 0;
	    double total_shallow_calls = 0;
	    AlleleFrequencyEstimate[] shallow_calls = null;
	    double EM_alt_freq = 0;
	    double EM_N = 0;
        double shallow_calls_fraction_correct = 0;

	    // 5. Call the pool.  (this step is the EM init)
	    pooled_caller.N = (int)true_N;
	    pooled_call = pooled_caller.map(tracker, ref, pooled_context);
	    System.out.print("POOLED_CALL " + pooled_call.asTabularString());

        // (this loop is the EM cycle)
	    EM_alt_freq = pooled_call.qstar; //pooled_call.qhat;
        int num_iterations = 10;
        double[] trajectory = new double[num_iterations + 1]; trajectory[0] = EM_alt_freq;
        double[] likelihood_trajectory = new double[num_iterations + 1]; likelihood_trajectory[0] = pooled_call.pBest;
        for (int iterations = 0; iterations < num_iterations; iterations++)
        { 
	        // 6. Re-call from shallow coverage using the estimated frequency as a prior, 
	        //    and compare to true deep calls, 
	        //    and compute new MAF estimate.
	        correct_shallow_calls = 0;
	        total_shallow_calls   = 0;
	        shallow_calls = new AlleleFrequencyEstimate[sample_names.size()];
	        EM_N        = 0.0;
            double EM_sum = 0.0;
            double likelihood = 0.0;

	        for (int i = 0; i < deep_calls.length; i++)
	        {
	            // Only shallow-call things where we know the truth!
	            if ((deep_calls[i].lodVsNextBest >= 5.0) || (deep_calls[i].lodVsRef <= -5.0))
	            {
	                shallow_callers.get(i).setAlleleFrequencyPrior(EM_alt_freq);
	                shallow_calls[i] = shallow_callers.get(i).map(tracker, ref, downsampled_contexts[i]);
	                String deep_genotype = deep_calls[i].genotype();
	                String shallow_genotype = shallow_calls[i].genotype();

                    likelihood += shallow_calls[i].lodVsNextBest;
	
	                //System.out.printf("DBG: %f %f %f %f\n", 
	                //                        deep_calls[i].lodVsNextBest,
	                //                        deep_calls[i].lodVsRef,
	                //                        shallow_calls[i].lodVsNextBest,
	                //                        shallow_calls[i].lodVsRef);

	                if (deep_genotype.equals(shallow_genotype)) 
                    { 
                        correct_shallow_calls += 1; 
                    }
	                total_shallow_calls += 1;
	                   
                    if (! FRACTIONAL_COUNTS)
                    {
			            EM_sum += shallow_calls[i].emperical_allele_frequency() * shallow_calls[i].N;
			            EM_N   += shallow_calls[i].N;
                    }
                    else
                    {
	                    for (int j = 0; j <= shallow_calls[i].N; j++)
	                    {
                            if (Double.isInfinite(shallow_calls[i].posteriors[j])) { shallow_calls[i].posteriors[j] = -10000; }
	                        System.out.printf("DBG3: %d %f %d\n", j, shallow_calls[i].posteriors[j], shallow_calls[i].N);
	                        EM_sum += Math.pow(10,shallow_calls[i].posteriors[j]) * (double)j;
	                        EM_N   += shallow_calls[i].N;
	                    }
                    }   

	            }
	        }
	        EM_alt_freq = EM_sum / EM_N;
	        shallow_calls_fraction_correct = correct_shallow_calls / total_shallow_calls;
            trajectory[iterations+1] = EM_alt_freq;
            likelihood_trajectory[iterations+1] = likelihood/(double)total_shallow_calls;

            System.out.printf("DBGTRAJ %f %f %f %f\n", EM_sum, EM_N, trajectory[iterations], trajectory[iterations+1]);
        }

        // 7. Compare to estimation from the pool.
        System.out.printf("EVAL %s %f %f %f %f %f %f %d %d %f %f %f %f\n",
                                pooled_call.location, 
                                pooled_call.lodVsRef, 
                                true_alt_freq, 
                                pooled_call.qhat, 
                                pooled_call.qstar, 
                                true_alt_freq * true_N,
                                pooled_call.emperical_allele_frequency() * true_N,
                                pooled_call.N, 
                                pooled_call.depth,
                                total_shallow_calls, 
                                correct_shallow_calls, 
                                shallow_calls_fraction_correct,
                                EM_alt_freq);
        for (int i = 0; i < likelihood_trajectory.length; i++)
        {
            System.out.printf("TRAJECTORY %f %f\n", trajectory[i], likelihood_trajectory[i]);
        }
        System.out.print("\n\n");

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
            String sample = read.getHeader().getReadGroup(RG).getSample();
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
