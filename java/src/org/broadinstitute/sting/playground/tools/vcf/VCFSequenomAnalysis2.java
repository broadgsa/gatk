/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.tools.vcf;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.*;
import java.util.*;
import java.lang.*;

import net.sf.picard.util.Interval;


class VCFSequenomAnalysis2 extends CommandLineProgram 
{
		@Argument(fullName = "sequenom", shortName = "sequenom", doc = "file to open", required = true) public String filename1;
		@Argument(fullName = "sequencing", shortName = "sequencing", doc = "file to open", required = true) public String filename2;
		@Argument(fullName = "out", shortName = "out", doc = "file to write results to", required = false) public String output_filename = "/dev/stdout";
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = true;
		@Argument(fullName = "verbose", shortName = "verbose", doc = "print extremely detailed stats", required = false) public Boolean verbose = false;
		@Argument(fullName = "qual_threshold", shortName = "qual_threshold", doc = "minimum genotype quality to consider", required = false) public long qual_threshold = 1;


		@Override
		protected int execute() 
		{
			//System.out.println("Loading " + filename + "...");
		
			PrintStream output = null;
			try
			{
				output = new PrintStream(new FileOutputStream(output_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			output.printf("PROBE interval flag ref alt n_total_sequenom failure_rate_sequenom n_alt_sequencing n_alt_sequenom p_alt_sequencing p_alt_sequenom HWE_sequencing_chi HWE_sequenom_chi HWE_sequencing_p HWE_sequenom_p is_singleton_in_sequencing singleton_matched_in_sequenom n_sequencing_hets n_sequenom_hets num_hets_also_het_in_sequenom num_hets_dropped_in_sequenom\n");

			VCFReader reader1;
			VCFReader reader2;

			if (autocorrect) 
			{ 
				reader1 = new VCFReader(VCFHomogenizer.create(filename1)); 
				reader2 = new VCFReader(VCFHomogenizer.create(filename2)); 
			}
			else 
			{ 
				reader1 = new VCFReader(new File(filename1)); 
				reader2 = new VCFReader(new File(filename2)); 
			}
		
			VCFHeader header1 = reader1.getHeader();
			VCFHeader header2 = reader2.getHeader();

			VCFRecord record1 = reader1.next();
			VCFRecord record2 = reader2.next();

			int[] sequenom_aaf_counts = new int[1024];
			int[] sequencing_aaf_counts = new int[1024];
			int max_aaf = 0;

			while(true)
			{
				if (record1 == null) { break; }
				if (record2 == null) { break; }

				Interval interval1 = VCFTool.getIntervalFromRecord(record1);
				Interval interval2 = VCFTool.getIntervalFromRecord(record2);

				int comparison = interval1.compareTo(interval2);
				
				if (comparison == 0)
				{
					// records match! compute concordance.
					
					// (unless one of them is "filtered")
					if (record1.isFiltered() || record2.isFiltered())
					{
						record1 = reader1.next();
						record2 = reader2.next();
						continue;
					}

					String[] sample_names = record2.getSampleNames();
					
					char ref = record1.getReference().charAt(0);
					char alt = VCFTool.getAlt(record2);
					
					int n_total_sequenom = VCFTool.Compute_n_total(record1, sample_names);
					int n_total_sequencing = VCFTool.Compute_n_total(record2, sample_names);
					double failure_rate_sequenom = VCFTool.Compute_failure_rate(record1);

					int n_alt_sequenom = VCFTool.Compute_n_alt(record1, sample_names);
					int n_alt_sequencing = VCFTool.Compute_n_alt(record2, sample_names);

					double p_alt_sequenom = (double)n_alt_sequenom / (double)n_total_sequenom;
					double p_alt_sequencing = (double)n_alt_sequencing / (double)n_total_sequencing;

					int n_het_sequenom = VCFTool.Compute_n_het(record1, sample_names);
					int n_het_sequencing = VCFTool.Compute_n_het(record2, sample_names);

					sequenom_aaf_counts[n_alt_sequenom] += 1;
					sequencing_aaf_counts[n_alt_sequencing] += 1;
					if (n_alt_sequenom > max_aaf) { max_aaf = n_alt_sequenom; }
					if (n_alt_sequencing > max_aaf) { max_aaf = n_alt_sequencing; }

					double HWE_sequenom = VCFTool.Compute_HWE(record1, sample_names);
					double HWE_sequencing = VCFTool.Compute_HWE(record2, sample_names);

					boolean isPolymorphic_sequenom   = (n_alt_sequenom   > 0) ? true : false;
					boolean isPolymorphic_sequencing = (n_alt_sequencing > 0) ? true : false;

					int is_singleton_in_sequencing = 0;
					int singleton_matched_in_sequenom = 0;
					
					if ((n_alt_sequencing == 1) && (VCFTool.Compute_n_alt(record1) > 0))
					{
						is_singleton_in_sequencing = 1;
						singleton_matched_in_sequenom = CheckSingletonMatch(record2, record1);
					}

					int[] het_match_ans = ComputeHetMatches(record2, record1);
					int num_hets_also_het_in_sequenom = het_match_ans[0];
				   	int num_hets_dropped_in_sequenom = het_match_ans[1];

					String flag = null;
					if (isPolymorphic_sequenom) { flag = "TP"; }
					else { flag = "FP"; }

					output.printf("PROBE %s %s %c %c %d %f %d %d %f %f %f %f %f %f %d %d %d %d %d %d\n",
										interval1,
										flag,
										ref,
										alt,
										n_total_sequenom,
										failure_rate_sequenom,
										n_alt_sequencing,
										n_alt_sequenom,
										p_alt_sequencing,
										p_alt_sequenom,
										HWE_sequencing,
										HWE_sequenom,
										VCFTool.P_from_Chi(HWE_sequencing),
										VCFTool.P_from_Chi(HWE_sequenom),
										is_singleton_in_sequencing,
										singleton_matched_in_sequenom,
										n_het_sequencing,
										n_het_sequenom,
										num_hets_also_het_in_sequenom,
				   						num_hets_dropped_in_sequenom);

					record1 = reader1.next();
					record2 = reader2.next();
				}
				else if (comparison > 0)
				{
					// interval1 is later than interval2.
					//System.err.printf("Skipping (2): %s\n", VCFTool.getIntervalFromRecord(record2));	
					record2 = reader2.next();
				}
				else if (comparison < 0)
				{
					// interval2 is later than interval1.
					//System.err.printf("Skipping (1): %s\n", VCFTool.getIntervalFromRecord(record1));	
					record1 = reader1.next();
				}

			}

			for (int i = 0; i < max_aaf; i++)
			{
				output.printf("AAF %d %d %d\n", i, sequenom_aaf_counts[i], sequencing_aaf_counts[i]);	
			}


			output.flush();
			output.close();


			return 0;
		}

		int CheckSingletonMatch(VCFRecord sequencing, VCFRecord sequenom)
		{
			String singleton_name = "";

			// first, check sequencing
			String[] sample_names = sequencing.getSampleNames();
			List<VCFGenotypeRecord> genotypes = sequencing.getVCFGenotypeRecords();
			int n_ref = 0;
			int n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
				if (g.charAt(1) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
			}
			if (n_alt != 1) { throw new RuntimeException(); }
			if (singleton_name.equals("")) { throw new RuntimeException(); }

			// now, check sequenom
			sample_names = sequenom.getSampleNames();
			genotypes = sequenom.getVCFGenotypeRecords();
			n_ref = 0;
			n_alt = 0;
			for (int i = 0; i < sample_names.length; i++)
			{
				if (sample_names[i].equals(singleton_name))
				{
					VCFGenotypeRecord rec = genotypes.get(i);
					List<VCFGenotypeEncoding> alleles = rec.getAlleles();
					String g = "";
					for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
					char[] c = g.toCharArray();
					Arrays.sort(c);
					g = new String(c);
					if (g.equals("..")) { continue; }
					if (g.charAt(0) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
					if (g.charAt(1) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; singleton_name = sample_names[i]; }
					break;
				}
			}
			if (n_alt > 0) { return 1; }
			else if (n_ref != 0) { return 0; }
			else { return -1; }
		}

		int[] ComputeHetMatches(VCFRecord sequencing, VCFRecord sequenom)
		{
			// first, check sequencing
			String[] sample_names = sequencing.getSampleNames();
			List<VCFGenotypeRecord> genotypes = sequencing.getVCFGenotypeRecords();
			ArrayList<String> het_samples = new ArrayList<String>();
			for (int i = 0; i < sample_names.length; i++)
			{
				int n_ref = 0;
				int n_alt = 0;

				VCFGenotypeRecord rec = genotypes.get(i);
				List<VCFGenotypeEncoding> alleles = rec.getAlleles();
				String g = "";
				for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
				char[] c = g.toCharArray();
				Arrays.sort(c);
				g = new String(c);
				if (g.equals("..")) { continue; }
				if (g.charAt(0) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
				if (g.charAt(1) == sequencing.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
				if (n_alt == 1) { het_samples.add(sample_names[i]); }
			}

			// now, check sequenom
			sample_names = sequenom.getSampleNames();
			genotypes = sequenom.getVCFGenotypeRecords();
			int matched_hets = 0;
			int dropped_hets = 0;
			int mismatched_hets = 0;
			int num_hets = het_samples.size();
			for (int i = 0; i < sample_names.length; i++)
			{
				if (het_samples.contains(sample_names[i]))
				{
					het_samples.remove(sample_names[i]);

					int n_ref = 0;
					int n_alt = 0;
					VCFGenotypeRecord rec = genotypes.get(i);
					List<VCFGenotypeEncoding> alleles = rec.getAlleles();
					String g = "";
					for (int j = 0; j < alleles.size(); j++) { g += alleles.get(j).getBases(); }
					char[] c = g.toCharArray();
					Arrays.sort(c);
					g = new String(c);
					if (g.equals("..")) { dropped_hets += 1; continue; }
					if (g.charAt(0) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
					if (g.charAt(1) == sequenom.getReference().charAt(0)) { n_ref += 1; } else { n_alt += 1; }
					if (n_alt == 1) { matched_hets += 1; }	
					else { mismatched_hets += 1; }
				}
			}

			if ((matched_hets + dropped_hets + mismatched_hets) != num_hets)
			{
				String warning = String.format("WARNING: %d + %d + %d != %d ", 
															matched_hets,
															dropped_hets,
															mismatched_hets,
															num_hets);
				for (int i = 0; i < het_samples.size(); i++) { warning += het_samples.get(i) + " "; }
				System.out.println(warning);
			}

			int[] ans = new int[2];
			ans[0] = matched_hets;
			ans[1] = dropped_hets;
			return ans;
		}
}
