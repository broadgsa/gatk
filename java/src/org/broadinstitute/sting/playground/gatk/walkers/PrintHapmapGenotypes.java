package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.gatk.walkers.varianteval.Histogram;

import java.util.*;
import java.io.*;

// Sanity check to test HapmapGenotypeROD
// Compute %dbsnp and transition/transversion rate.

@By(DataSource.REFERENCE)
@Requires(DataSource.REFERENCE)
@Allows(DataSource.REFERENCE)
public class PrintHapmapGenotypes extends RefWalker<Integer, Integer> 
{
    //@Argument(required=false, shortName="n_frequency_bins", doc="") public int n_frequency_bins = 20;

    public void initialize() 
	{
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) 
	{
        // Iterate over each analysis, and update it
        //rodDbSNP dbsnp =   (rodDbSNP)tracker.lookup("dbsnp", null);
        HapMapGenotypeROD A = (HapMapGenotypeROD)tracker.lookup("A", null);

		if (A != null) 
		{ 
			GenomeLoc loc = A.getLocation();
			String[] sample_ids = A.getSampleIDs();
			String[] genotypes  = A.getGenotypes();

			for (int i = 0; i < sample_ids.length; i++)
			{
				out.printf("%s %s %s\n", loc, sample_ids[i], genotypes[i]);
			}
			out.printf("\n");
		}	

		return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) 
	{
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) 
	{
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) 
	{
    }

}
