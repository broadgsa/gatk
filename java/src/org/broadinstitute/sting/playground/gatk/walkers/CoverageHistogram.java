
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;

import java.util.*;

// Plot a histogram of depth of coverage
// j.maguire 6-11-2009

@By(DataSource.REFERENCE)
public class CoverageHistogram extends LocusWalker<Integer,Integer>
{
	//@Argument(fullName="start", shortName="start", required=false, doc="start") public Integer START = 0;

	// Private state.
	long[] coverage_hist;
	int max_depth;

	long sum_coverage;
	long num_sites;

	/////////
	// Walker Interface Functions 
    public void initialize() 
	{
		coverage_hist = new long[1000000];
        Arrays.fill(coverage_hist, 0);
		max_depth = 0;

		sum_coverage = 0;
		num_sites = 0;
	}

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
	{
		if (ref.getBase() == 'N') { return null; }
		int depth = context.getReads().size();
		coverage_hist[depth] += 1;
		if (depth > max_depth) { max_depth = depth; }

		sum_coverage += depth;
		num_sites += 1;

		return 0;
	}

    public void onTraversalDone(Integer sum) 
	{
        double mean_coverage = (double)sum_coverage / (double)num_sites;
        double mean_good_coverage = (double)sum_coverage / ((double)(num_sites - coverage_hist[0]));
		out.printf("# all_sites                  : mean:%f num_sites:%d%n", mean_coverage, num_sites);
        out.printf("# sites with at least 1 read : mean:%f num_sites:%d%n", mean_good_coverage, num_sites - coverage_hist[0]);

        // Code for calculting std devs adapted from Michael Melgar's python script

        // Find the maximum extent of 'good' data
        // First, find the mode
        long maxValue = coverage_hist[1]; // ignore doc=0
        int mode = 1;
        for (int i = 2; i <= max_depth; i++) {
            if (coverage_hist[i] > maxValue) {
                maxValue = coverage_hist[i];
                mode = i;
            }
        }
        // now, procede to find end of good Gaussian fit
        long dist = (long)Math.pow(10, 9);
        while ( Math.abs(coverage_hist[mode] - coverage_hist[1]) < dist )
            dist = Math.abs(coverage_hist[mode++] - coverage_hist[1]);
        int maxGoodDepth = mode + 1;

        // calculate the mean of the good region
        long totalGoodSites = 0, totalGoodDepth = 0;
        for (int i = 1; i <= maxGoodDepth; i++) { // ignore doc=0
            totalGoodSites += coverage_hist[i];
            totalGoodDepth += i * coverage_hist[i];
        }
        double meanGoodDepth = (double)totalGoodDepth / (double)totalGoodSites;

        // calculate the variance and standard deviation of the good region
        double var = 0.0;
        for (int i = 1; i <= maxGoodDepth; i++) {  // ignore doc=0
            var += coverage_hist[i] * Math.pow(meanGoodDepth - (double)i, 2);
        }
        double stdev = Math.sqrt(var / (double)totalGoodSites);
        out.printf("# sites within Gaussian fit  : mean:%f num_sites:%d std_dev:%f%n", meanGoodDepth, totalGoodSites, stdev);

        for (int i = 1; i <= 5; i++)
            out.printf("# Gaussian mean + %d Std Dev  : %f%n", i, (meanGoodDepth + i*stdev));

		out.println("\ndepth count freq(percent)");
		for (int i = 0; i <= max_depth; i++)
		{
			out.printf("%d %d %f\n", i, coverage_hist[i], (100.0*coverage_hist[i]) / (double)num_sites);
		}
		return;
	}

    public Integer reduceInit() 
	{
		return 0;
	}

    public Integer reduce(Integer record, Integer sum) 
	{
		return 0;
	}

	// END Walker Interface Functions
	/////////
}
