
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.util.zip.*;
import java.io.*;

// Beta iterative multi-sample caller
// j.maguire 6-11-2009

public class MultiSampleCallerAccuracyTest extends MultiSampleCaller
{
    @Argument(required=false, shortName="lod_threshold", doc="") public double LOD_THRESHOLD = 1e-6;


    public void initialize() 
	{
		this.DISCOVERY_OUTPUT = "/dev/null";
		this.INDIVIDUAL_OUTPUT = "/dev/null";
		super.initialize();
	}

    public MultiSampleCallResult map(RefMetaDataTracker tracker, char ref, LocusContext context) 
	{
        HapMapGenotypeROD hapmap = (HapMapGenotypeROD)tracker.lookup("hapmap", null);

		MultiSampleCallResult call_result = super.map(tracker, ref, context);
		EM_Result em_result = call_result.em_result;

		// Compute individual accuracy.
		double n_calls = 0;
		double n_correct = 0;
		for (int i = 0; i < em_result.sample_names.length; i++)
		{
			String sample_name = em_result.sample_names[i];
			String hyp_genotype = em_result.genotype_likelihoods[i].BestGenotype();
			String ref_genotype = hapmap.get(sample_name);
			double lod = em_result.genotype_likelihoods[i].LodVsNextBest();

			if ((lod > LOD_THRESHOLD) && (ref_genotype != null))
			{
				n_calls += 1;
				if (hyp_genotype.equals(ref_genotype))
				{
					n_correct += 1;
				}
			}

		}

		out.printf("%s %.0f %.0f %.2f%%\n", 
						context.getLocation(), 
						n_calls, 
						n_correct, 
						100.0*n_correct / n_calls);

		return call_result;
	}

    public void onTraversalDone(String sum) 
	{
		out.println("MultiSampleCallerAccuracyTest done.");
		return;
	}

    public String reduceInit() 
	{
		return super.reduceInit();
	}

    public String reduce(MultiSampleCallResult record, String sum) 
	{
		return super.reduce(record, sum);
	}

	// END Walker Interface Functions
	/////////

}
