/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ÓSoftwareÓ), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ÓAS ISÓ, WITHOUT WARRANTY OF ANY KIND,
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

import net.sf.picard.util.Interval;


class VCFSequenomAnalysis extends CommandLineProgram 
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

			output.printf("interval flag ref alt missing_base n_total_sequenom failure_rate_sequenom n_alt_sequencing HWE_sequencing_chi HWE_sequenom_chi HWE_sequencing_p HWE_sequenom_p\n");

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


			while(true)
			{
				if (record1 == null) { break; }
				if (record2 == null) { break; }

				String[] sample_names = record2.getSampleNames();

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
					
					char ref = record1.getReference().charAt(0);
					char alt = VCFTool.getAlt(record2);
					
					int n_total_sequenom = VCFTool.Compute_n_total(record1);
					double failure_rate_sequenom = VCFTool.Compute_failure_rate(record1);

					int n_alt_sequenom = VCFTool.Compute_n_alt(record1);
					int n_alt_sequencing = VCFTool.Compute_n_alt(record2);

					double HWE_sequenom = VCFTool.Compute_HWE(record1, sample_names);
					double HWE_sequencing = VCFTool.Compute_HWE(record2);

					boolean isPolymorphic_sequenom   = (n_alt_sequenom   > 0) ? true : false;
					boolean isPolymorphic_sequencing = (n_alt_sequencing > 0) ? true : false;

					String flag = null;
					char missing_base = '.';

					if (isPolymorphic_sequenom) 
					{ 
						flag = "TP"; 
						if (n_alt_sequenom == n_total_sequenom) { missing_base = ref; }
					}
					else 
					{ 
						flag = "FP"; 
						missing_base = alt;
					}

					output.printf("%s %s %c %c %c %d %f %d %f %f %f %f\n",
										interval1,
										flag,
										ref,
										alt,
										missing_base,
										n_total_sequenom,
										failure_rate_sequenom,
										n_alt_sequencing,
										HWE_sequencing,
										HWE_sequenom,
										VCFTool.P_from_Chi(HWE_sequencing),
										VCFTool.P_from_Chi(HWE_sequenom));

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

			output.flush();
			output.close();


			return 0;
		}
}
