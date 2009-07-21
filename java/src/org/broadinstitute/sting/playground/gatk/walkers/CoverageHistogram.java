
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.util.zip.*;
import java.io.*;

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

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) 
	{
		if (ref == 'N') { return null; }
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
		out.printf("# all_sites                  : mean:%f num_sites:%d%n", mean_coverage, num_sites);
        out.printf("# sites with at least 1 read : mean:%f num_sites:%d%n%n", sum_coverage / ((double)(num_sites - coverage_hist[0])), num_sites - coverage_hist[0]);

		out.println("depth count freq");
		for (int i = 0; i <= max_depth; i++)
		{
			out.printf("%d %d %f\n", i, coverage_hist[i], (100.0*coverage_hist[i]) / num_sites);
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
