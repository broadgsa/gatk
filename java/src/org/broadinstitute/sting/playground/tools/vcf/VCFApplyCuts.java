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

//import org.broadinstitute.sting.playground.tools.vcf.VCFOptimize.Cut;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;

import org.broadinstitute.sting.utils.genotype.vcf.*;


import java.io.*;
import java.util.*;

//import org.apache.commons.math.optimization.*;
//import org.apache.commons.math.optimization.direct.*;
//import org.apache.commons.math.analysis.MultivariateRealFunction;

// First draft of a program for working with VCF files in various ways.


/**
 * @author jmaguire
 */


class VCFApplyCuts extends CommandLineProgram 
{
		@Argument(fullName = "vcf", shortName = "vcf", doc = "file to open", required = true) public String filename;
		@Argument(fullName = "auto_correct", shortName = "auto_correct", doc = "auto-correct the VCF file if it's off-spec", required = false) public Boolean autocorrect = false;
		@Argument(fullName = "cuts", shortName = "cuts", doc = "file to read cuts from", required = true) public String cuts_filename;
		@Argument(fullName = "output", shortName = "output", doc = "file to write filtered VCF to", required = true) public String output_filename;

			class Cut
			{
				public double lod;
				public double slod;
				public int freq;

				public Cut(double lod, double slod)
				{
					this.lod = lod;
					this.slod = slod;
					this.freq = -1;
				}

				public Cut(double lod, double slod, int freq)
				{
					this.lod = lod;
					this.slod = slod;
					this.freq = freq;
				}

				public Cut(String record)
				{
					String[] tokens = record.split("\\s+");
					this.freq = Integer.parseInt(tokens[0]);					
					this.lod  = Double.parseDouble(tokens[1]);					
					this.slod = Double.parseDouble(tokens[2]);					
				}

				public String toString()
				{
					return String.format("%d %f %f", freq, lod, slod);
				}
			}

			private boolean applyCuts(ArrayList<Cut> cuts, VCFHeader header, VCFRecord record)
			{
				Map<String,String> info = record.getInfoValues();
				
				if (! info.containsKey("AC")) 
				{ 
					throw new RuntimeException("AC not present in record: \n" + record.toStringEncoding(header));
				}
				if (! info.containsKey("DP")) 
				{ 
					throw new RuntimeException("DP not present in record: \n" + record.toStringEncoding(header));
				}
				if (! info.containsKey("SB")) 
				{ 
					throw new RuntimeException("SB not present in record: \n" + record.toStringEncoding(header));
				}
	
				boolean transition = VCFTool.isTransition(record);
				int freq           = Integer.parseInt(record.getInfoValues().get("AC"));
				double LOD         = record.getQual();
				double depth       = Double.parseDouble(record.getInfoValues().get("DP"));
				double SLOD        = Double.parseDouble(record.getInfoValues().get("SB"));

				for (int i = 0; i < cuts.size(); i++)
				{
					Cut cut = cuts.get(i);
					if (cut.freq == freq)
					{
						if ((LOD >= cut.lod) && (-1*SLOD >= cut.slod)) { return true; }
					}
				}

				return false;
			}


		@Override
		protected int execute() 
		{
			// load cuts.

			ArrayList<Cut> cuts = null;
			Scanner cuts_file = null;
			try
			{
				cuts = new ArrayList<Cut>();
				cuts_file = new Scanner(new File(cuts_filename));
			}
			catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			while (cuts_file.hasNextLine())
			{
				String line = cuts_file.nextLine();
				Cut cut = new Cut(line);
				cuts.add(cut);
			}

		
			VCFReader reader = null;

			if (autocorrect) { reader = new VCFReader(new File(filename),new VCFHomogenizer()); }
			else { reader = new VCFReader(new File(filename)); }

			VCFHeader header = reader.getHeader();


			VCFWriter writer = new VCFWriter(new File(output_filename));
            writer.writeHeader(header);


			Date start_time = new Date();
			int n_records_processed = 0;
			int n_records_passed = 0;
			while(reader.hasNext())
			{
				VCFRecord record = reader.next();

				if (applyCuts(cuts, header, record) == true)
				{
					writer.addRecord(record);
					n_records_passed += 1;
				}	
		
				n_records_processed += 1;
			}
			System.out.printf("Processed %d records\n", n_records_processed);
			System.out.printf("Passed %d records\n", n_records_passed);

			writer.close();
			
			return 0;
		}
}
