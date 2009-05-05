
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisTK;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class PoolCaller extends LocusWalker<AlleleFrequencyEstimate, String> 
{
    List<SingleSampleGenotyper> callers = null;
    List<String> sample_names = null;

    //@Argument(required=false, shortName="log_metrics", defaultValue="true") public boolean LOG_METRICS;
    @Argument(required=false, shortName="fractional_counts", doc="fractional counts", defaultValue="false") public boolean FRACTIONAL_COUNTS;
    
    private Random random;

    private SAMFileHeader header;

    public void initialize() 
    { 
        GenomeAnalysisTK toolkit = this.getToolkit();
        this.header = toolkit.getEngine().getSAMHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();

        /*
        GenomeAnalysisTK toolkit = this.getToolkit();
        SAMFileHeader header = toolkit.getSamReader().getFileHeader();
        List<SAMReadGroupRecord> read_groups = header.getReadGroups();
        */

        sample_names    = new ArrayList<String>();
        callers         = new ArrayList<SingleSampleGenotyper>();

        random = new Random(42);

        for (int i = 0; i < read_groups.size(); i++)
        {
            String sample_name = read_groups.get(i).getSample();
            sample_names.add(sample_name);
            //System.out.println("SAMPLE: " + sample_name);

            SingleSampleGenotyper caller = new SingleSampleGenotyper();
            caller.metricsFileName = "/dev/null";
            caller.lodThreshold = 5.0;
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
            //System.out.printf("DEPTH %s %d\n", sample_names.get(i), contexts[i].getReads().size());
        }
        //System.out.printf("DEPTH %s %d\n", "TOTAL", context.getReads().size());

        // EM Loop:
	    AlleleFrequencyEstimate[] calls = null;
	    double EM_alt_freq = 0;
	    double EM_N = 0;

        // (this loop is the EM cycle)
	    EM_alt_freq = 0.5;
        int num_iterations = 10;
        double[] trajectory = new double[num_iterations + 1]; trajectory[0] = EM_alt_freq;
        double[] likelihood_trajectory = new double[num_iterations + 1]; likelihood_trajectory[0] = 0.0;
        for (int iterations = 0; iterations < num_iterations; iterations++)
        { 
	        // 6. Re-call from shallow coverage using the estimated frequency as a prior, 
	        //    and compare to true deep calls, 
	        //    and compute new MAF estimate.
	        calls = new AlleleFrequencyEstimate[sample_names.size()];
	        EM_N        = 0.0;
            double EM_sum = 0.0;
            double likelihood = 0.0;

	        for (int i = 0; i < sample_names.size(); i++)
	        {
	                callers.get(i).setAlleleFrequencyPrior(EM_alt_freq);
	                calls[i] = callers.get(i).map(tracker, ref, contexts[i]);
	                String genotype = calls[i].genotype();

                    likelihood += calls[i].posterior();
	
	                //System.out.printf("DBG: %f %f %f %f\n", 
	                //                        deep_calls[i].lodVsNextBest,
	                //                        deep_calls[i].lodVsRef,
	                //                        shallow_calls[i].lodVsNextBest,
	                //                        shallow_calls[i].lodVsRef);
	
                    if (! FRACTIONAL_COUNTS)
                    {
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

        for (int i = 0; i < sample_names.size(); i++)
        {
            ReadBackedPileup pileup = new ReadBackedPileup(ref, contexts[i]);  
            System.out.printf("CALL %s %s %f %f %s\n", sample_names.get(i), calls[i].genotype(), calls[i].lodVsRef, calls[i].lodVsNextBest, pileup.getBases());
        }

        // 7. Compare to estimation from the pool.
        System.out.printf("EVAL %s %f\n",
                                context.getLocation(),
                                EM_alt_freq);
        //for (int i = 0; i < likelihood_trajectory.length; i++)
        //{
        //    System.out.printf("TRAJECTORY %f %f\n", trajectory[i], likelihood_trajectory[i]);
        //}
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
