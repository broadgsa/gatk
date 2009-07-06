
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
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

public class CoverageHistogram extends LocusWalker<Integer,Integer>
{
	//@Argument(fullName="start", shortName="start", required=false, doc="start") public Integer START = 0;

	// Private state.
	int[] coverage_hist;
	int max_depth;

	long sum_coverage;
	long num_sites;

	/////////
	// Walker Interface Functions 
    public void initialize() 
	{
		coverage_hist = new int[1000000];
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
		out.printf("# mean:%f num_sites:%d\n\n", mean_coverage, num_sites);

		out.println("depth count");
		for (int i = 1; i < max_depth; i++)
		{
			out.printf("%d %d\n", i, coverage_hist[i]);
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
