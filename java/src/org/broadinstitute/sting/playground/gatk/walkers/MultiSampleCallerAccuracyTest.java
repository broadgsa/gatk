
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
import org.broadinstitute.sting.playground.indels.Matrix;

import java.util.*;
import java.util.zip.*;
import java.io.*;

// Beta iterative multi-sample caller
// j.maguire 6-11-2009

public class MultiSampleCallerAccuracyTest extends MultiSampleCaller
{
    @Argument(required=false, shortName="lod_threshold", doc="") public double LOD_THRESHOLD = 1e-6;
    @Argument(required=true, shortName="stats_output", doc="") public String STATS_OUTPUT;

	Matrix<Integer> n_variants;
	Matrix<Integer> n_found;

	PrintStream stats_output;

    public void initialize() 
	{
		this.DISCOVERY_OUTPUT = "/dev/null";
		this.INDIVIDUAL_OUTPUT = "/dev/null";

		super.initialize();

		n_variants = new Matrix<Integer>(sample_names.size()*2, sample_names.size()*2);
		n_found    = new Matrix<Integer>(sample_names.size()*2, sample_names.size()*2);

		for (int i = 0; i < sample_names.size()*2; i++)
		{
			for (int j = 0; j < sample_names.size()*2; j++)
			{
				n_variants.set(i,j,0);
				n_found.set(i,j,0);
			}
		}

		try
		{
			stats_output = new PrintStream(STATS_OUTPUT);
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}

	}

    public MultiSampleCallResult map(RefMetaDataTracker tracker, char ref, LocusContext context) 
	{
        HapMapGenotypeROD hapmap = (HapMapGenotypeROD)tracker.lookup("hapmap", null);

		// Collect all the variants and the normals.
		ArrayList<String> variant_samples = new ArrayList<String>();
		ArrayList<String> reference_samples = new ArrayList<String>();

		int n_ref_chromosomes = 0;
		int n_alt_chromosomes = 0;

		String reference_genotype = String.format("%c%c", ref, ref);
		for (int i = 0; i < sample_names.size(); i++)
		{
			String true_genotype = hapmap.get(sample_names.get(i));
			if (true_genotype == null) { continue; }

			if (true_genotype.equals(reference_genotype)) { reference_samples.add(sample_names.get(i)); }
			else { variant_samples.add(sample_names.get(i)); }

			if (true_genotype.equals(reference_genotype)) { n_ref_chromosomes += 1; }
			else if (true_genotype.contains(String.format("%c",ref))) { n_ref_chromosomes += 1; n_alt_chromosomes += 1; }
			else { n_alt_chromosomes += 2; }

		}

			// Put together a context.
			ArrayList<String> working_samples = new ArrayList<String>();
			working_samples.addAll(variant_samples);
			working_samples.addAll(reference_samples);
			LocusContext working_context = filterLocusContextBySamples(context, working_samples);

			// Call.
			MultiSampleCallResult call_result = super.map(tracker, ref, working_context);
			EM_Result em_result = call_result.em_result;

			// Compute Statistics.
			if (n_variants == null) { System.out.printf("n_variants is null\n"); }
			if (n_found == null) { System.out.printf("n_found is null\n"); }
			n_variants.set(n_ref_chromosomes, n_alt_chromosomes, n_variants.get(n_ref_chromosomes, n_alt_chromosomes)+1);
			if ((call_result.lod > LOD_THRESHOLD) && (n_alt_chromosomes >= 1))
			{
				n_found.set(n_ref_chromosomes, n_alt_chromosomes, n_found.get(n_ref_chromosomes, n_alt_chromosomes)+1);
			}

		return null;
	}

	private void PrintStats()
	{
		stats_output.printf("n_reference_chromosomes n_variant_chromosomes n_sites n_found fraction_found\n");
		for (int i = 0; i < sample_names.size()*2; i++)
		{
			for (int j = 0; j < sample_names.size()*2; j++)
			{
				int N = (int)n_variants.get(i,j);
				int found = (int)n_found.get(i,j);

				if (N == 0) { continue; }
				if (found == 0) { continue; }

				double fraction_found = 100.0 * (double)found / (double)N;
				n_variants.set(i,j,0);
				n_found.set(i,j,0);
				stats_output.printf("%d %d %d %d %f\n", 
										i,
										j,
										N,
										found,
										fraction_found);
			}
		}
	}

    public void onTraversalDone(String sum) 
	{
		PrintStats();
		stats_output.flush();
		stats_output.close();
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


	/////////
	// BEGIN Utility Functions

	// Filter a locus context by sample IDs
	//   (pulls out only reads from the specified samples, and returns them in one context).
    private LocusContext filterLocusContextBySamples(LocusContext context, List<String> sample_names)
    {
		HashSet<String> index = new HashSet<String>();
		for (int i = 0; i < sample_names.size(); i++)
		{
			index.add(sample_names.get(i));
		}

		ArrayList<SAMRecord> reads = new ArrayList();
		ArrayList<Integer> offsets = new ArrayList();

        for (int i = 0; i < context.getReads().size(); i++)
        {
            SAMRecord read = context.getReads().get(i);
            Integer offset = context.getOffsets().get(i);
            String RG = (String)(read.getAttribute("RG"));

            assert(header != null);
            assert(header.getReadGroup(RG) != null);

            String sample = header.getReadGroup(RG).getSample();
			if (SAMPLE_NAME_REGEX != null) { sample = sample.replaceAll(SAMPLE_NAME_REGEX, "$1"); }

			if (index.contains(sample))
			{
            	reads.add(read); 
            	offsets.add(offset); 
			}
        }

		return new LocusContext(context.getLocation(), reads, offsets);
    }

	// END Utility Functions
	/////////

}
